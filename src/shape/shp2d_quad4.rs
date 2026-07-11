use crate::shape::shape_base::*;

pub struct Quad4;

impl Shape2D for Quad4 {
    fn num_quad(&self) -> usize {
        4
    }
    fn w(&self) -> Vec<f64> {
        vec![1.0, 1.0, 1.0, 1.0]
    }
    fn a(&self) -> Vec<f64> {
        vec![-0.5773502691896257, -0.5773502691896257, 0.5773502691896257, 0.5773502691896257]
    }
    fn b(&self) -> Vec<f64> {
        vec![-0.5773502691896257, 0.5773502691896257, -0.5773502691896257, 0.5773502691896257]
    }

    fn num_node(&self) -> usize {
        4
    }
    fn eval(&self, a: f64, b: f64) -> Vec<f64> {
        vec![
            0.25 * (1.0 - a) * (1.0 - b),
            0.25 * (1.0 + a) * (1.0 - b),
            0.25 * (1.0 + a) * (1.0 + b),
            0.25 * (1.0 - a) * (1.0 + b)
        ]
    }
    fn grad(&self, a: f64, b: f64) -> (Vec<f64>, Vec<f64>) {
        let dn_da = vec![
            -0.25 * (1.0 - b),
            0.25 * (1.0 - b),
            0.25 * (1.0 + b),
            -0.25 * (1.0 + b)
        ];
        let dn_db = vec![
            -0.25 * (1.0 - a),
            -0.25 * (1.0 + a),
            0.25 * (1.0 + a), 
            0.25 * (1.0 - a)
        ];
        (dn_da, dn_db)
    }
}
