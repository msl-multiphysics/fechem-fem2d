// lin2 constants
pub const W_LIN2: [f64; 2] = [1.0, 1.0];
pub const A_LIN2: [f64; 2] = [-0.5773502691896257, 0.5773502691896257];

// lin2 shape function
pub fn lin2_eval(a: f64) -> [f64; 2] {
    [0.5 * (1.0 - a), 0.5 * (1.0 + a)]
}

// lin2 shape function gradients
pub fn lin2_grad(_a: f64) -> [f64; 2] {
    [-0.5, 0.5]
}
