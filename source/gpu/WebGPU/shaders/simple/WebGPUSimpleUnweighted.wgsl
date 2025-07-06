const inv_width: f32 = 4;

alias atomic_f32 = atomic<u32>;
alias Bin = atomic_f32;

struct Atom {
    xyz: vec3<f32>,
    w: f32,
}

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

fn atomic_add_Bin(bin: ptr<storage,Bin,read_write>, value: f32) {
    atomic_add_f32(bin, value);
}

@group(0) @binding(0) var<storage,read> atom_buffer_1: array<Atom>;
@group(0) @binding(1) var<storage,read> atom_buffer_2: array<Atom>;
@group(0) @binding(2) var<storage,read_write> histogram: array<Bin>;

@compute @workgroup_size(64)
fn calculate_self(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms = arrayLength(&atom_buffer_1);
    if (id.x >= num_atoms) {
        return;
    }

    // consider using a local histogram
    // var local_histogram: array<atomic<u32>, arrayLength(&histogram)>;
    let atom1 = atom_buffer_1[id.x];
    atomic_add_Bin(&histogram[0], atom1.w*atom1.w);
    for (var i: u32 = id.x + 1u; i < num_atoms; i = i + 1u) {
        let atom2 = atom_buffer_1[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w * atom2.w;
        atomic_add_Bin(&histogram[bin], weight);
    }
}

@compute @workgroup_size(64)
fn calculate_cross(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms1 = arrayLength(&atom_buffer_1);
    let num_atoms2 = arrayLength(&atom_buffer_2); 
    if (id.x >= num_atoms1) {
        return;
    }

    let atom1 = atom_buffer_1[id.x];
    for (var i: u32 = 0; i < num_atoms2; i = i + 1u) {
        let atom2 = atom_buffer_2[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w*atom2.w;
        atomic_add_Bin(&histogram[bin], weight);
    }
}