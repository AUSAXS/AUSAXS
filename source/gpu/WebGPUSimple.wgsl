const var inv_width = 0.1;

@group(0) @binding(0) var<storage,read> inputBuffer: array<f32,64>;
@group(0) @binding(1) var<storage,read_write> outputBuffer: array<f32,64>;

@compute @workgroup_size(32)
fn eval_single(array<f32, 4> v1, array<f32, 4> v2) -> array<f32, 2> {
    return array<f32, 2>(round(inv_width*distance(v1, v2)), v1[3]*v2[3])
}