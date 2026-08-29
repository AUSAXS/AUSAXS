// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// Unweighted distance histogram. Mirrors hist::distance_calculator::SimpleCPU with weighted_bins = false.
//
// Two things make a naive implementation unusable here. Accumulating directly into the global histogram
// costs a contended device-memory atomic per pair distance, which dominates the runtime by an order of
// magnitude, so contributions are first gathered in a workgroup-local histogram and flushed once per
// workgroup. And accumulating in f32 is far too imprecise: the bins reach values around 1e11, where the
// single-precision spacing is several thousand and adding a term of order 1 does nothing at all. The
// bins are therefore fixed-point integers, spread over two u32 words that are added with the wrapping
// atomicAdd of the low word carrying into the high word. That makes the sum exact and independent of
// the order the atomics happen to be applied in, so results are reproducible between runs.
//
// The local histogram only covers the first LOCAL_BINS bins, as that is what fits in workgroup storage.
// Longer distances are rare in practice and are added directly to the global histogram.

const WORKGROUP_SIZE: u32 = 256u;
const LOCAL_BINS: u32 = 2048u; // 16 kB, the guaranteed workgroup storage size

// Fixed-point scale of the bin values. Must match ausaxs::gpu::fixed_point_scale.
const FIXED_SCALE: f32 = 65536.0;
const TWO_32: f32 = 4294967296.0;

struct Atom {
    xyz: vec3<f32>,
    w: f32,
}

// A bin value as a 64-bit fixed-point integer. Must match ausaxs::gpu::Buffers<false>::HistogramType.
struct Bin {
    low: atomic<u32>,
    high: atomic<u32>,
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
@group(0) @binding(2) var<storage,read_write> histogram: array<Bin>;
@group(0) @binding(3) var<uniform> params: Params;

var<workgroup> local_low: array<atomic<u32>, LOCAL_BINS>;
var<workgroup> local_high: array<atomic<u32>, LOCAL_BINS>;

// Split a value into the low and high word of its signed 64-bit fixed-point representation. The
// subtraction is exact, so the pair represents the scaled value to the full precision of the input.
fn to_fixed(value: f32) -> vec2<u32> {
    let scaled = round(value*FIXED_SCALE);
    let high = floor(scaled/TWO_32);
    return vec2<u32>(u32(scaled - high*TWO_32), bitcast<u32>(i32(high)));
}

fn add_local(bin: u32, value: vec2<u32>) {
    let old = atomicAdd(&local_low[bin], value.x);
    var high = value.y;
    if (old + value.x < old) { // the low word wrapped, so carry into the high word
        high = high + 1u;
    }
    if (high != 0u) {
        atomicAdd(&local_high[bin], high);
    }
}

fn add_global(bin: u32, value: vec2<u32>) {
    let old = atomicAdd(&histogram[bin].low, value.x);
    var high = value.y;
    if (old + value.x < old) {
        high = high + 1u;
    }
    if (high != 0u) {
        atomicAdd(&histogram[bin].high, high);
    }
}

// distances beyond the axis range are discarded rather than folded into the last bin
fn add_to_bin(bin: u32, weight: f32) {
    if (bin >= params.bin_count) {
        return;
    }
    if (bin < LOCAL_BINS) {
        add_local(bin, to_fixed(weight));
    } else {
        add_global(bin, to_fixed(weight));
    }
}

fn clear_local(local_index: u32) {
    for (var bin = local_index; bin < LOCAL_BINS; bin = bin + WORKGROUP_SIZE) {
        atomicStore(&local_low[bin], 0u);
        atomicStore(&local_high[bin], 0u);
    }
}

fn flush_local(local_index: u32) {
    let bins = min(LOCAL_BINS, params.bin_count);
    for (var bin = local_index; bin < bins; bin = bin + WORKGROUP_SIZE) {
        let value = vec2<u32>(atomicLoad(&local_low[bin]), atomicLoad(&local_high[bin]));
        if (value.x != 0u || value.y != 0u) {
            add_global(bin, value);
        }
    }
}

// all distances from atom i to the atoms after it
fn evaluate_self_row(i: u32, num_atoms: u32, scale: f32) {
    let atom1 = atom_buffer_1[i];
    add_to_bin(0u, scale*atom1.w*atom1.w); // the self-distance skipped by the loop below

    // the pair (i, j) is only visited once, so every distance is counted twice
    let w1 = 2*scale*atom1.w;
    for (var j: u32 = i + 1u; j < num_atoms; j = j + 1u) {
        let atom2 = atom_buffer_1[j];
        let d = distance(atom1.xyz, atom2.xyz);
        add_to_bin(u32(round(params.inv_width*d)), w1*atom2.w);
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
            add_to_bin(u32(round(params.inv_width*d)), w1*atom2.w);
        }
    }

    workgroupBarrier();
    flush_local(local_index);
}
