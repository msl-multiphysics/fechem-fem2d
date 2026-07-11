use crate::shape::shape_base::*;

pub struct Tri3;

impl Shape2D for Tri3 {
    fn num_quad(&self) -> usize {
        3
    }
    fn w(&self) -> Vec<f64> {
        vec![1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]
    }
    fn a(&self) -> Vec<f64> {
        vec![1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]
    }
    fn b(&self) -> Vec<f64> {
        vec![1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0]
    }

    fn num_node(&self) -> usize {
        3
    }
    fn eval(&self, a: f64, b: f64) -> Vec<f64> {
        vec![
            1.0 - a - b,
            a,
            b
        ]
    }
    fn grad(&self, _a: f64, _b: f64) -> (Vec<f64>, Vec<f64>) {
        let dn_da = vec![
            -1.0,
            1.0,
            0.0
        ];
        let dn_db = vec![
            -1.0,
            0.0,
            1.0
        ];
        (dn_da, dn_db)
    }
}
