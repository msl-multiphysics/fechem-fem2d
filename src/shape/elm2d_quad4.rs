// quad4 constants
pub const M_FRAC_1_SQRT3: f64 = 0.5773502691896257; // 1/sqrt(3)
pub const W_QUAD4: [f64; 4] = [1.0, 1.0, 1.0, 1.0];
pub const A_QUAD4: [f64; 4] = [-M_FRAC_1_SQRT3, -M_FRAC_1_SQRT3, M_FRAC_1_SQRT3, M_FRAC_1_SQRT3];
pub const B_QUAD4: [f64; 4] = [-M_FRAC_1_SQRT3, M_FRAC_1_SQRT3, -M_FRAC_1_SQRT3, M_FRAC_1_SQRT3];

// quad4 shape function
pub fn quad4_eval(a: f64, b: f64) -> [f64; 4] {
    [
        0.25 * (1.0 - a) * (1.0 - b),
        0.25 * (1.0 + a) * (1.0 - b),
        0.25 * (1.0 + a) * (1.0 + b),
        0.25 * (1.0 - a) * (1.0 + b),
    ]
}

// quad4 shape function gradients
pub fn quad4_grad(a: f64, b: f64) -> ([f64; 4], [f64; 4]) {
    let dn_da = [
        -0.25 * (1.0 - b),
        0.25 * (1.0 - b),
        0.25 * (1.0 + b),
        -0.25 * (1.0 + b),
    ];
    let dn_db = [
        -0.25 * (1.0 - a),
        -0.25 * (1.0 + a),
        0.25 * (1.0 + a),
        0.25 * (1.0 - a),
    ];
    (dn_da, dn_db)
}
