use crate::shape::shape_base::*;

pub struct Lin2;

impl Shape1D for Lin2 {
    fn num_quad(&self) -> usize {
        2
    }
    fn w(&self) -> Vec<f64> {
        vec![1.0, 1.0]
    }
    fn a(&self) -> Vec<f64> {
        vec![-0.5773502691896257, 0.5773502691896257]
    }

    fn num_node(&self) -> usize {
        2
    }
    fn eval(&self, a: f64) -> Vec<f64> {
        vec![
            0.5 * (1.0 - a),
            0.5 * (1.0 + a)
        ]
    }
    fn grad(&self, a: f64) -> Vec<f64> {
        vec![
            -0.5,
            0.5
        ]
    }
}
