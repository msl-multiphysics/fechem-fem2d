// tri3 constants
pub const W_TRI3: [f64; 1] = [0.5];
pub const A_TRI3: [f64; 1] = [1.0 / 3.0];
pub const B_TRI3: [f64; 1] = [1.0 / 3.0];

// tri3 shape function
pub fn tri3_eval(a: f64, b: f64) -> [f64; 3] {
    [1.0 - a - b, a, b]
}

// tri3 shape function gradients
pub fn tri3_grad(_a: f64, _b: f64) -> ([f64; 3], [f64; 3]) {
    let dn_da = [-1.0, 1.0, 0.0];
    let dn_db = [-1.0, 0.0, 1.0];
    (dn_da, dn_db)
}
