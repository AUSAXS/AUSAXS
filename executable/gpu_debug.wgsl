const inv_width: f32 = 4;

struct Atom {
    xyz: vec3<f32>,
    w: f32,
};

@group(0) @binding(0) var<storage,read> atom_buffer_1: array<Atom>;
@group(0) @binding(1) var<storage,read> atom_buffer_2: array<Atom>;
@group(0) @binding(2) var<storage,read_write> histogram: array<atomic<u32>>;

fn f(x: f32) -> f32 {
    return 2.0 * x + 1.0;
}

// workaround until wgsl supports atomic operations on f32
fn atomic_add(i: u32, value: f32) {
    var old = atomicLoad(&histogram[i]);
    loop {
        let new_value = value + bitcast<f32>(old);
        let exchange_result = atomicCompareExchangeWeak(&histogram[i], old, bitcast<u32>(new_value));
        if exchange_result.exchanged {
            return;
        }
        old = exchange_result.old_value;
    }
}

@compute @workgroup_size(32)
fn calculate_self(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms = arrayLength(&atom_buffer_1);

    if (id.x >= num_atoms) {
        return;
    }

    let atom1 = atom_buffer_1[id.x];
    atomic_add(0, atom1.w*atom1.w);
    for (var i = id.x + 1u; i < num_atoms; i = i + 1u) {
        let atom2 = atom_buffer_1[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w * atom2.w;
        atomic_add(bin, weight);
    }
}

@compute @workgroup_size(32)
fn calculate_cross(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms1 = arrayLength(&atom_buffer_1);
    let num_atoms2 = arrayLength(&atom_buffer_2); 

    if (id.x >= num_atoms1) {
        return;
    }

    let atom1 = atom_buffer_1[id.x];
    for (var i = num_atoms1; i < num_atoms1 + num_atoms2; i = i + 1u) {
        let atom2 = atom_buffer_2[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w*atom2.w;
        atomic_add(bin, weight);
    }
}

@compute @workgroup_size(2)
fn calculate_test(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms = arrayLength(&atom_buffer_1);

    if (id.x >= num_atoms) {
        return;
    }

    let atom1 = atom_buffer_1[id.x];
    atomicAdd(&histogram[id.x], u32(atom1.w));
}