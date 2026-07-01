use crate::base::geom_bnd::Boundary;
use crate::base::error::FEChemError;

#[derive(Clone, Copy, PartialEq)]
pub enum ScalarBoundaryType {
    Constant,
    Function,
}

pub struct ScalarBoundary {
    // ids
    pub sclbnd_id: usize,
    pub bnd_id: usize,  // boundary this scalar is attached to

    // values
    pub node_value: Vec<f64>,  // [nid] -> values at nodes
    pub quad_value: Vec<Vec<f64>>, // [eid][qid] -> values at quadrature points
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied

    // scalar data
    pub scl_type: ScalarBoundaryType,
    pub con_value: f64,
    pub fun_func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64>,  // time, [x, y], [scalars]
    pub fun_sclid: Vec<usize>,  // ids of scalars used in function
}

impl Default for ScalarBoundary{
    fn default() -> ScalarBoundary {
        ScalarBoundary {
            sclbnd_id: 0,
            bnd_id: 0,
            node_value: Vec::new(),
            quad_value: Vec::new(),
            node_dir: Vec::new(),
            scl_type: ScalarBoundaryType::Constant,
            con_value: 0.0,
            fun_func: Box::new(|_: f64, _: [f64; 2], _: &[f64]| 0.0),
            fun_sclid: vec![],
        }
    }
}

impl ScalarBoundary {
    pub fn new_from_constant(sclbnd_id: usize, bnd: &Boundary, value_const: f64) -> Result<ScalarBoundary, FEChemError> {
        // create struct
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // set values
        sclbnd.node_value = vec![value_const; bnd.num_node];
        sclbnd.quad_value = Vec::new();
        for eid in 0..bnd.num_elem {
            sclbnd.quad_value.push(vec![value_const; bnd.elem_node[eid]]);
        }
        sclbnd.node_dir = vec![false; bnd.num_node];

        // set constant type
        sclbnd.scl_type = ScalarBoundaryType::Constant;
        sclbnd.con_value = value_const;

        // result
        Ok(sclbnd)
    }
}
