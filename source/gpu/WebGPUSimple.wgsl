const var inv_width = 4;

struct Atom {
    xyz: vec3<f32>;
    w: f32;
};

@group(0) @binding(0) var<storage,read> atom_buffer_1: array<Atom>;
@group(0) @binding(1) var<storage,read> atom_buffer_2: array<Atom>;
@group(0) @binding(2) var<storage,read_write> histogram: array<f32>;

@compute @workgroup_size(32)
fn calculate_self() {
    let global_id = @builtin(global_invocation_id).x;
    let num_atoms = arrayLength(&atom_buffer_1);

    if (global_id >= num_atoms) {
        return;
    }

    let atom1 = atom_buffer_1[global_id];

    for (var i = global_id + 1u; i < num_atoms; i = i + 1u) {
        let atom2 = atom_buffer_1[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = atom1.w * atom2.w;

        atomicAdd(&histogram[bin], weight);
    }
}

@compute @workgroup_size(32)
fn calculate_cross() {
    let global_id = @builtin(global_invocation_id).x;
    let num_atoms1 = arrayLength(&atom_buffer_1);
    let num_atoms2 = arrayLength(&atom_buffer_2); 

    if (global_id >= num_atoms1) {
        return;
    }

    let atom1 = atom_buffer_1[global_id];

    for (var i = num_atoms1; i < num_atoms1 + num_atoms2; i = i + 1u) {
        let atom2 = atom_buffer_2[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = atom1.w * atom2.w;

        atomicAdd(&histogram[bin], weight);
    }
}