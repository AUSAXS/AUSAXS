// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// Weighted distance histogram. Mirrors hist::distance_calculator::SimpleCPU with weighted_bins = true.
// In addition to the bin value, the number of contributions and the sum of the raw distances are
// tracked so the weighted bin center can be recovered on the host.

alias atomic_f32 = atomic<u32>;

struct Atom {
    xyz: vec3<f32>,
    w: f32,
}

// Must match ausaxs::gpu::Buffers<true>::HistogramType.
struct WeightedBin {
    value: atomic_f32,
    count: atomic<u32>,
    center: atomic_f32,
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

// workaround until wgsl supports atomic operations on f32
fn atomic_add_f32(v1: ptr<storage,atomic_f32,read_write>, v2: f32) {
    var old_v = atomicLoad(v1);
    loop {
        let new_v = v2 + bitcast<f32>(old_v);
        let exchange_result = atomicCompareExchangeWeak(v1, old_v, bitcast<u32>(new_v));
        if exchange_result.exchanged {
            return;
        }
        old_v = exchange_result.old_value;
    }
}

// distances beyond the axis range are discarded rather than folded into the last bin
fn add_to_bin(bin: u32, count: u32, weight: f32, distance: f32) {
    if (bin >= params.bin_count) {
        return;
    }
    atomicAdd(&histogram[bin].count, count);
    atomic_add_f32(&histogram[bin].value, weight);
    atomic_add_f32(&histogram[bin].center, distance);
}

@compute @workgroup_size(64)
fn calculate_self(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms = arrayLength(&atom_buffer_1);
    if (id.x >= num_atoms) {
        return;
    }
    let scale = f32(params.scale);

    // the self-distance skipped by the loop below. the CPU kernel tracks its weight as the bin count,
    // which is harmless as the bin center of a zero distance is zero regardless.
    let self_weight = scale*atom_buffer_1[id.x].w*atom_buffer_1[id.x].w;
    add_to_bin(0u, u32(self_weight), self_weight, 0.0);

    // the pair (i, j) is only visited once, so every distance is counted twice
    let atom1 = atom_buffer_1[id.x];
    let w1 = 2*scale*atom1.w;
    for (var i: u32 = id.x + 1u; i < num_atoms; i = i + 1u) {
        let atom2 = atom_buffer_1[i];
        let d = distance(atom1.xyz, atom2.xyz);
        add_to_bin(u32(round(params.inv_width*d)), 2*params.scale, w1*atom2.w, 2*scale*d);
    }
}

@compute @workgroup_size(64)
fn calculate_cross(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms1 = arrayLength(&atom_buffer_1);
    let num_atoms2 = arrayLength(&atom_buffer_2);
    if (id.x >= num_atoms1) {
        return;
    }
    let scale = f32(params.scale);

    let atom1 = atom_buffer_1[id.x];
    let w1 = 2*scale*atom1.w;
    for (var i: u32 = 0; i < num_atoms2; i = i + 1u) {
        let atom2 = atom_buffer_2[i];
        let d = distance(atom1.xyz, atom2.xyz);
        add_to_bin(u32(round(params.inv_width*d)), 2*params.scale, w1*atom2.w, 2*scale*d);
    }
}
