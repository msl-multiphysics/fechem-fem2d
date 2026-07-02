use crate::base::geom_bnd::Boundary;
use crate::base::error::FEChemError;
use crate::base::itg_bnd::IntegralBoundary;
use crate::base::vars::Variables;
use crate::base::write_csv::write_sclbnd_csv;
use crate::base::write_vtu::write_sclbnd_vtu;
use crate::shape::prelude::*;
use faer::Col;
use faer::linalg::solvers::Solve;
use faer::sparse::Triplet;
use faer::sparse::SparseColMat;
use faer::sparse::linalg::solvers::{Lu, SymbolicLu};

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
    // constant - fixes node and quad values
    // function - updates quad values per iteration; updates node values only when writing
    pub node_value: Vec<f64>,  // [nid] -> values at nodes
    pub quad_value: Vec<Vec<f64>>, // [eid][qid] -> values at quadrature points
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied

    // scalar data
    pub scl_type: ScalarBoundaryType,
    pub con_value: f64,
    pub fun_func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64>,  // time, [x, y], [scalars]
    pub fun_sclid: Vec<usize>,  // ids of scalars used in function

    // output file
    pub file_name: String,  // path to file without extension
    pub file_type: String,  // leave empty if no output
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
            file_name: String::new(),
            file_type: String::new(),
        }
    }
}

impl ScalarBoundary {
    pub fn new_from_constant(sclbnd_id: usize, bnd: &Boundary, value_const: f64, file_path: String) -> Result<ScalarBoundary, FEChemError> {
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

        // set outputs if file path is not empty
        if file_path == "" {
            sclbnd.file_name = String::new();
            sclbnd.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            sclbnd.file_name = parts[0..parts.len()-1].join(".");
            sclbnd.file_type = parts[parts.len()-1].to_string();
        }

        // result
        Ok(sclbnd)
    }

    pub fn new_from_function(sclbnd_id: usize, bnd: &Boundary, value_func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64>, scldom_ids: Vec<usize>, file_path: String) -> Result<ScalarBoundary, FEChemError> {
        // create struct
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // set values
        // actual values will be computed later
        sclbnd.node_value = vec![0.0; bnd.num_node];
        sclbnd.quad_value = Vec::new();
        for eid in 0..bnd.num_elem {
            sclbnd.quad_value.push(vec![0.0; bnd.elem_node[eid]]);
        }
        sclbnd.node_dir = vec![false; bnd.num_node];

        // set function type
        sclbnd.scl_type = ScalarBoundaryType::Function;
        sclbnd.fun_func = value_func;
        sclbnd.fun_sclid = scldom_ids;

        // set outputs if file path is not empty
        if file_path == "" {
            sclbnd.file_name = String::new();
            sclbnd.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            sclbnd.file_name = parts[0..parts.len()-1].join(".");
            sclbnd.file_type = parts[parts.len()-1].to_string();
        }

        // result
        Ok(sclbnd)
    }

    pub fn write(&self, bnd: &Boundary, ts: usize) -> Result<(), FEChemError> {
        // // write depending on file type
        match self.file_type.as_str() {
            "csv" => write_sclbnd_csv(&bnd, &self, ts)?,
            "vtu" => write_sclbnd_vtu(&bnd, &self, ts)?,
            "" => (),  // do nothing if file type is empty
            _ => panic!("Unsupported file type: {}", self.file_type),
        }

        // placeholder
        Ok(())
    }

    pub fn transfer_quad_node(&mut self, bnd: &Boundary, itgbnd: &IntegralBoundary) {
        // skip if not function type
        if self.scl_type != ScalarBoundaryType::Function {
            return;
        }
        
        // initialize global matrix and vector
        let num_node = bnd.num_node;
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let mut b_vec = Col::zeros(num_node);

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // initialize local matrices
            let num_node_loc = bnd.elem_node[eid];
            let mut m_loc = vec![vec![0.0; num_node_loc]; num_node_loc];
            let mut b_loc = vec![0.0; num_node_loc];

            // get integral data
            let num_quad = itgbnd.num_quad[eid];
            let jac_det = &itgbnd.jac_det[eid];

            // get properties
            let quad_val = &self.quad_value[eid];

            // assemble local mass matrix and rhs vector
            match num_node_loc {
                2 => {
                    for qid in 0..num_quad {
                        let n = lin2_eval(A_LIN2[qid]);
                        let coeff = W_LIN2[qid] * jac_det[qid].sqrt();
                        for v in 0..num_node_loc {
                            for j in 0..num_node_loc {
                                m_loc[v][j] += coeff * n[v] * n[j];
                            }
                            b_loc[v] += coeff * quad_val[qid] * n[v];
                        }
                    }
                }
                _ => {panic!("Invalid element type");}
            }

            // assemble global matrix and vector

            // iterate over local matrix entries
            let node_id = &bnd.elem_node_id[eid];
            for v in 0..num_node_loc {
                for j in 0..num_node_loc {
                    // get node ids
                    let nid_v = node_id[v];
                    let nid_j = node_id[j];

                    // add to global matrix
                    a_triplet.push(Triplet::new(nid_v, nid_j, m_loc[v][j]));
                }
            }
            for v in 0..num_node_loc {
                // get node ids
                let nid_v = node_id[v];

                // add to global vector
                b_vec[nid_v] += b_loc[v];
            }

        }

        // solve mass system
        let m_mat = SparseColMat::try_new_from_triplets(num_node, num_node, &a_triplet)
            .map_err(|_| FEChemError::FailedMatrixSolve {caller: "ScalarBoundary::transfer_quad_node".to_string()}).unwrap();
        let symbolic = SymbolicLu::try_new(m_mat.symbolic())
            .map_err(|_| FEChemError::FailedMatrixSolve {caller: "ScalarBoundary::transfer_quad_node".to_string()}).unwrap();
        let lu = Lu::try_new_with_symbolic(symbolic, m_mat.as_ref())
            .map_err(|_| FEChemError::FailedMatrixSolve {caller: "ScalarBoundary::transfer_quad_node".to_string()}).unwrap();
        let u_vec = lu.solve(&b_vec);

        // store nodal values
        for nid in 0..num_node {
            self.node_value[nid] = u_vec[nid];
        }

    }
}

// this was separated from the struct to prevent borrowing issues
pub fn sclbnd_update_function(vars: &mut Variables, sclbnd_id: usize, t: f64) {
    if vars.scl_bnd[sclbnd_id].scl_type != ScalarBoundaryType::Function {
        return;
    }

    // get objects
    let bnd_id = vars.scl_bnd[sclbnd_id].bnd_id;
    let bnd = &vars.bnd[bnd_id];
    let itgbnd = &vars.itg_bnd[bnd_id];

    // iterate over quadrature points
    for eid in 0..bnd.num_elem {
        for qid in 0..itgbnd.num_quad[eid] {
            let value = sclbnd_quad_value(vars, sclbnd_id, eid, qid, t);
            vars.scl_bnd[sclbnd_id].quad_value[eid][qid] = value;
        }
    }
}

pub fn sclbnd_quad_value(vars: &Variables, sclbnd_id: usize, eid: usize, qid: usize, t: f64) -> f64 {
    let sclbnd = &vars.scl_bnd[sclbnd_id];
    match sclbnd.scl_type {
        ScalarBoundaryType::Constant => sclbnd.con_value,
        ScalarBoundaryType::Function => {
            let itgbnd = &vars.itg_bnd[sclbnd.bnd_id];
            let x = itgbnd.quad_x[eid][qid];
            let y = itgbnd.quad_y[eid][qid];
            let mut scl_val: Vec<f64> = Vec::new();
            for &sclbnd_sub in &sclbnd.fun_sclid {
                scl_val.push(sclbnd_quad_value(vars, sclbnd_sub, eid, qid, t));
            }
            (sclbnd.fun_func)(t, [x, y], scl_val.as_slice())
        }
    }
}

