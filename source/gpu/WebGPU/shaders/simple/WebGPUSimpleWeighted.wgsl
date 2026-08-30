// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// Weighted distance histogram. Mirrors hist::distance_calculator::SimpleCPU with weighted_bins = true.
// In addition to the bin value, the number of contributions and the sum of the raw distances are
// tracked so the weighted bin center can be recovered on the host.
//
// See WebGPUSimpleUnweighted.wgsl for why the accumulation goes through a workgroup-local histogram of
// fixed-point integers. Each bin needs three counters here, so fewer of them fit in workgroup storage.

const WORKGROUP_SIZE: u32 = 256u;
const LOCAL_BINS: u32 = 800u; // 16 kB, the guaranteed workgroup storage size

// Fixed-point scale of the bin values. Must match ausaxs::gpu::fixed_point_scale.
const FIXED_SCALE: f32 = 65536.0;
const TWO_32: f32 = 4294967296.0;

struct Atom {
    xyz: vec3<f32>,
    w: f32,
}

// Must match ausaxs::gpu::Buffers<true>::HistogramType. The value and center are 64-bit fixed-point
// integers spread over two u32 words each.
struct WeightedBin {
    value_low: atomic<u32>,
    value_high: atomic<u32>,
    count: atomic<u32>,
    center_low: atomic<u32>,
    center_high: atomic<u32>,
}

// Per-dispatch parameters. Must match ausaxs::gpu::ShaderParameters.
struct Params {
    inv_width: f32,  // inverse bin width; bin = round(inv_width*distance)
    scale: u32,      // multiplicative factor applied to every contribution
    bin_count: u32,  // number of bins in the histogram buffer
    padding: u32,
}

@group(0) @binding(0) var<storage,read> atom_buffer_1: array<Atom>;
@group(0) @binding(1) var<storage,read> atom_buffer_2: array<Atom>;
@group(0) @binding(2) var<storage,read_write> histogram: array<WeightedBin>;
@group(0) @binding(3) var<uniform> params: Params;

var<workgroup> local_value_low: array<atomic<u32>, LOCAL_BINS>;
var<workgroup> local_value_high: array<atomic<u32>, LOCAL_BINS>;
var<workgroup> local_count: array<atomic<u32>, LOCAL_BINS>;
var<workgroup> local_center_low: array<atomic<u32>, LOCAL_BINS>;
var<workgroup> local_center_high: array<atomic<u32>, LOCAL_BINS>;

// Split a value into the low and high word of its signed 64-bit fixed-point representation. The
// subtraction is exact, so the pair represents the scaled value to the full precision of the input.
fn to_fixed(value: f32) -> vec2<u32> {
    let scaled = round(value*FIXED_SCALE);
    let high = floor(scaled/TWO_32);
    return vec2<u32>(u32(scaled - high*TWO_32), bitcast<u32>(i32(high)));
}

// wgsl rejects passing two pointers into the same variable to one function, so the low and high words
// cannot be forwarded as arguments and each accumulator needs its own copy of the carrying add.
fn add_local_value(bin: u32, value: vec2<u32>) {
    let old = atomicAdd(&local_value_low[bin], value.x);
    var carried = value.y;
    if (old + value.x < old) { // the low word wrapped, so carry into the high word
        carried = carried + 1u;
    }
    if (carried != 0u) {
        atomicAdd(&local_value_high[bin], carried);
    }
}

fn add_local_center(bin: u32, value: vec2<u32>) {
    let old = atomicAdd(&local_center_low[bin], value.x);
    var carried = value.y;
    if (old + value.x < old) {
        carried = carried + 1u;
    }
    if (carried != 0u) {
        atomicAdd(&local_center_high[bin], carried);
    }
}

fn add_global_value(bin: u32, value: vec2<u32>) {
    let old = atomicAdd(&histogram[bin].value_low, value.x);
    var carried = value.y;
    if (old + value.x < old) {
        carried = carried + 1u;
    }
    if (carried != 0u) {
        atomicAdd(&histogram[bin].value_high, carried);
    }
}

fn add_global_center(bin: u32, value: vec2<u32>) {
    let old = atomicAdd(&histogram[bin].center_low, value.x);
    var carried = value.y;
    if (old + value.x < old) {
        carried = carried + 1u;
    }
    if (carried != 0u) {
        atomicAdd(&histogram[bin].center_high, carried);
    }
}

// distances beyond the axis range are discarded rather than folded into the last bin
fn add_to_bin(bin: u32, count: u32, weight: f32, center: f32) {
    if (bin >= params.bin_count) {
        return;
    }
    if (bin < LOCAL_BINS) {
        atomicAdd(&local_count[bin], count);
        add_local_value(bin, to_fixed(weight));
        add_local_center(bin, to_fixed(center));
    } else {
        atomicAdd(&histogram[bin].count, count);
        add_global_value(bin, to_fixed(weight));
        add_global_center(bin, to_fixed(center));
    }
}

fn clear_local(local_index: u32) {
    for (var bin = local_index; bin < LOCAL_BINS; bin = bin + WORKGROUP_SIZE) {
        atomicStore(&local_value_low[bin], 0u);
        atomicStore(&local_value_high[bin], 0u);
        atomicStore(&local_count[bin], 0u);
        atomicStore(&local_center_low[bin], 0u);
        atomicStore(&local_center_high[bin], 0u);
    }
}

fn flush_local(local_index: u32) {
    let bins = min(LOCAL_BINS, params.bin_count);
    for (var bin = local_index; bin < bins; bin = bin + WORKGROUP_SIZE) {
        let count = atomicLoad(&local_count[bin]);
        if (count == 0u) {
            continue;
        }
        atomicAdd(&histogram[bin].count, count);
        add_global_value(bin, vec2<u32>(atomicLoad(&local_value_low[bin]), atomicLoad(&local_value_high[bin])));
        add_global_center(bin, vec2<u32>(atomicLoad(&local_center_low[bin]), atomicLoad(&local_center_high[bin])));
    }
}

// all distances from atom i to the atoms after it
fn evaluate_self_row(i: u32, num_atoms: u32, scale: f32) {
    let atom1 = atom_buffer_1[i];

    // the self-distance skipped by the loop below. the CPU kernel tracks its weight as the bin count,
    // which is harmless as the bin center of a zero distance is zero regardless.
    let self_weight = scale*atom1.w*atom1.w;
    add_to_bin(0u, u32(self_weight), self_weight, 0.0);

    // the pair (i, j) is only visited once, so every distance is counted twice
    let w1 = 2*scale*atom1.w;
    for (var j: u32 = i + 1u; j < num_atoms; j = j + 1u) {
        let atom2 = atom_buffer_1[j];
        let d = distance(atom1.xyz, atom2.xyz);
        add_to_bin(u32(round(params.inv_width*d)), 2*params.scale, w1*atom2.w, 2*scale*d);
    }
}

@compute @workgroup_size(WORKGROUP_SIZE)
fn calculate_self(@builtin(global_invocation_id) id: vec3<u32>, @builtin(local_invocation_index) local_index: u32) {
    clear_local(local_index);
    workgroupBarrier();

    // row i has num_atoms-i-1 distances to evaluate, so each invocation takes one row from either end
    // to give all of them the same amount of work. note that the barriers above and below must be
    // reached by every invocation, so no early return is possible here.
    let num_atoms = arrayLength(&atom_buffer_1);
    if (num_atoms != 0u && id.x < (num_atoms + 1u)/2u) {
        let scale = f32(params.scale);
        evaluate_self_row(id.x, num_atoms, scale);

        let mirrored = num_atoms - 1u - id.x;
        if (mirrored != id.x) {
            evaluate_self_row(mirrored, num_atoms, scale);
        }
    }

    workgroupBarrier();
    flush_local(local_index);
}

@compute @workgroup_size(WORKGROUP_SIZE)
fn calculate_cross(@builtin(global_invocation_id) id: vec3<u32>, @builtin(local_invocation_index) local_index: u32) {
    clear_local(local_index);
    workgroupBarrier();

    let num_atoms1 = arrayLength(&atom_buffer_1);
    let num_atoms2 = arrayLength(&atom_buffer_2);
    if (id.x < num_atoms1) {
        let scale = f32(params.scale);
        let atom1 = atom_buffer_1[id.x];
        let w1 = 2*scale*atom1.w;
        for (var j: u32 = 0; j < num_atoms2; j = j + 1u) {
            let atom2 = atom_buffer_2[j];
            let d = distance(atom1.xyz, atom2.xyz);
            add_to_bin(u32(round(params.inv_width*d)), 2*params.scale, w1*atom2.w, 2*scale*d);
        }
    }

    workgroupBarrier();
    flush_local(local_index);
}
