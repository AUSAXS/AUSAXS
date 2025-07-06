const inv_width: f32 = 4;

struct Atom {
    xyz: vec3<f32>,
    w: f32,
}

struct atomic_f32 {
    value: atomic<u32>
}

struct Bin {
    value: atomic_f32
}

struct WeightedBin {
    value: atomic_f32,
    count: atomic<u32>,
    center: atomic_f32
}

// workaround until wgsl supports atomic operations on f32
fn atomic_add_f32(v1: ptr<storage,atomic_f32,read_write>, v2: f32) {
    var old_v = atomicLoad(&v1.value);
    loop {
        let new_v = v2 + bitcast<f32>(old_v);
        let exchange_result = atomicCompareExchangeWeak(&v1.value, old_v, bitcast<u32>(new_v));
        if exchange_result.exchanged {
            return;
        }
        old_v = exchange_result.old_value;
    }
}

fn atomic_add_Bin(bin: ptr<storage,Bin,read_write>, value: f32) {
    atomic_add_f32(&bin.value, value);
}

fn atomic_add_WeightedBin(bin: ptr<storage,WeightedBin,read_write>, value: f32, distance: f32) {
    atomicAdd(&bin.count, 1);
    atomic_add_f32(&bin.value, value);
    atomic_add_f32(&bin.center, distance);
}

//###########################################//
//###          Unweighted shaders         ###//
//###########################################//
@group(0) @binding(0) var<storage,read> unweighted_atom_buffer_1: array<Atom>;
@group(0) @binding(1) var<storage,read> unweighted_atom_buffer_2: array<Atom>;
@group(0) @binding(2) var<storage,read_write> unweighted_histogram: array<Bin>;

@compute @workgroup_size(64)
fn calculate_self(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms = arrayLength(&unweighted_atom_buffer_1);
    if (id.x >= num_atoms) {
        return;
    }

    // consider using a local histogram
    // var local_histogram: array<atomic<u32>, arrayLength(&histogram)>;
    let atom1 = unweighted_atom_buffer_1[id.x];
    atomic_add_Bin(&unweighted_histogram[0], atom1.w*atom1.w);
    for (var i: u32 = id.x + 1u; i < num_atoms; i = i + 1u) {
        let atom2 = unweighted_atom_buffer_1[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w * atom2.w;
        atomic_add_Bin(&unweighted_histogram[bin], weight);
    }
}

@compute @workgroup_size(64)
fn calculate_cross(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms1 = arrayLength(&unweighted_atom_buffer_1);
    let num_atoms2 = arrayLength(&unweighted_atom_buffer_2); 
    if (id.x >= num_atoms1) {
        return;
    }

    let atom1 = unweighted_atom_buffer_1[id.x];
    for (var i: u32 = 0; i < num_atoms2; i = i + 1u) {
        let atom2 = unweighted_atom_buffer_2[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w*atom2.w;
        atomic_add_Bin(&unweighted_histogram[bin], weight);
    }
}

//###########################################//
//###           Weighted shaders          ###//
//###########################################//
@group(1) @binding(0) var<storage,read> weighted_atom_buffer_1: array<Atom>;
@group(1) @binding(1) var<storage,read> weighted_atom_buffer_2: array<Atom>;
@group(1) @binding(2) var<storage,read_write> weighted_histogram: array<WeightedBin>;

@compute @workgroup_size(64)
fn calculate_weighted_self(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms = arrayLength(&weighted_atom_buffer_1);
    if (id.x >= num_atoms) {
        return;
    }

    let atom1 = weighted_atom_buffer_1[id.x];
    atomic_add_WeightedBin(&weighted_histogram[0], atom1.w*atom1.w, 0.0);
    for (var i: u32 = id.x + 1u; i < num_atoms; i = i + 1u) {
        let atom2 = weighted_atom_buffer_1[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w * atom2.w;
        atomic_add_WeightedBin(&weighted_histogram[bin], weight, distance);
    }
}

@compute @workgroup_size(64)
fn calculate_weighted_cross(@builtin(global_invocation_id) id: vec3<u32>) {
    let num_atoms1 = arrayLength(&weighted_atom_buffer_1);
    let num_atoms2 = arrayLength(&weighted_atom_buffer_2); 
    if (id.x >= num_atoms1) {
        return;
    }

    let atom1 = weighted_atom_buffer_1[id.x];
    for (var i: u32 = 0; i < num_atoms2; i = i + 1u) {
        let atom2 = weighted_atom_buffer_2[i];
        let distance = distance(atom1.xyz, atom2.xyz);
        let bin = u32(round(inv_width * distance));
        let weight = 2*atom1.w*atom2.w;
        atomic_add_WeightedBin(&weighted_histogram[bin], weight, distance);
    }
}