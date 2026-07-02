use crate::base::error::FEChemError;
use crate::base::geom_dom::Domain;
use crate::base::itg_dom::IntegralDomain;
use crate::base::write_csv::write_scldom_csv;
use crate::base::write_vtu::write_scldom_vtu;
use crate::base::vars::Variables;
use crate::shape::prelude::*;
use faer::Col;
use faer::linalg::solvers::Solve;
use faer::sparse::Triplet;
use faer::sparse::SparseColMat;
use faer::sparse::linalg::solvers::{Lu, SymbolicLu};

#[derive(Clone, Copy, PartialEq)]
pub enum ScalarDomainType {
    Constant,
    Function,
    Unknown,
}

pub struct ScalarDomain {
    // ids
    pub scldom_id: usize,
    pub dom_id: usize,  // domain this scalar is attached to

    // values
    // constant - fixes node and quad values
    // function - updates quad values per iteration; updates node values only when writing
    // unknown - updates node values per iteration; never updates quad values
    pub node_value: Vec<f64>,  // [nid] -> values at nodes
    pub quad_value: Vec<Vec<f64>>, // [eid][qid] -> values at quadrature points
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied

    // scalar type
    pub scl_type: ScalarDomainType,
    pub con_value: f64,
    pub fun_func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64>,  // time, [x, y], [scalars]
    pub fun_sclid: Vec<usize>,  // ids of scalars used in function
    pub unk_start: usize,  // index of this scalar in solution vector

    // output file
    pub file_name: String,  // path to file without extension
    pub file_type: String,  // leave empty if no output
}

impl Default for ScalarDomain{
    fn default() -> ScalarDomain {
        ScalarDomain {
            scldom_id: 0,
            dom_id: 0,
            node_value: Vec::new(),
            quad_value: Vec::new(),
            node_dir: Vec::new(),
            scl_type: ScalarDomainType::Constant,
            con_value: 0.0,
            fun_func: Box::new(|_: f64, _: [f64; 2], _: &[f64]| 0.0),
            fun_sclid: vec![],
            unk_start: 0,
            file_name: String::new(),
            file_type: String::new(),
        }
    }
}

impl ScalarDomain {
    pub fn new_from_constant(scldom_id: usize, dom: &Domain, value_const: f64, file_path: String) -> Result<ScalarDomain, FEChemError> {
        // create struct
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // set values
        scldom.node_value = vec![value_const; dom.num_node];
        scldom.quad_value = Vec::new();
        for eid in 0..dom.num_elem {
            scldom.quad_value.push(vec![value_const; dom.elem_node[eid]]);
        }
        scldom.node_dir = vec![false; dom.num_node];

        // set constant type
        scldom.scl_type = ScalarDomainType::Constant;
        scldom.con_value = value_const;

        // set outputs if file path is not empty
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            scldom.file_name = parts[0..parts.len()-1].join(".");
            scldom.file_type = parts[parts.len()-1].to_string();
        }

        // result
        Ok(scldom)
    }

    pub fn new_from_function(scldom_id: usize, dom: &Domain, value_func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64>, scldom_ids: Vec<usize>, file_path: String) -> Result<ScalarDomain, FEChemError> {
        // create struct
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // set values
        // actual values will be computed later
        scldom.node_value = vec![0.0; dom.num_node];
        scldom.quad_value = Vec::new();
        for eid in 0..dom.num_elem {
            scldom.quad_value.push(vec![0.0; dom.elem_node[eid]]);
        }
        scldom.node_dir = vec![false; dom.num_node];

        // set function type
        scldom.scl_type = ScalarDomainType::Function;
        scldom.fun_func = value_func;
        scldom.fun_sclid = scldom_ids;

        // set outputs if file path is not empty
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            scldom.file_name = parts[0..parts.len()-1].join(".");
            scldom.file_type = parts[parts.len()-1].to_string();
        }

        // result
        Ok(scldom)
    }

    pub fn new_from_unknown(scldom_id: usize, dom: &Domain, value_init: f64, file_path: String) -> Result<ScalarDomain, FEChemError> {
        // create struct
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // set values
        scldom.node_value = vec![value_init; dom.num_node];
        scldom.quad_value = Vec::new();
        for eid in 0..dom.num_elem {
            scldom.quad_value.push(vec![value_init; dom.elem_node[eid]]);
        }
        scldom.node_dir = vec![false; dom.num_node];

        // set unknown type
        scldom.scl_type = ScalarDomainType::Unknown;
        scldom.unk_start = 0;

        // set outputs if file path is not empty
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            scldom.file_name = parts[0..parts.len()-1].join(".");
            scldom.file_type = parts[parts.len()-1].to_string();
        }

        // result
        Ok(scldom)
    }

    pub fn write(&self, dom: &Domain, ts: usize) -> Result<(), FEChemError> {
        // write depending on file type
        match self.file_type.as_str() {
            "csv" => write_scldom_csv(&dom, &self, ts)?,
            "vtu" => write_scldom_vtu(&dom, &self, ts)?,
            "" => (),  // do nothing if file type is empty
            _ => panic!("Unsupported file type: {}", self.file_type),
        }

        // placeholder
        Ok(())
    }

    pub fn update_unknown(&mut self, dom: &Domain, x_vec: &Col<f64>) {
        // skip if not unknown type
        if self.scl_type != ScalarDomainType::Unknown {
            return;
        }

        // update nodal values
        let num_node = dom.num_node;
        for nid in 0..num_node {
            let xid = self.unk_start + nid;
            let value = x_vec[xid];
            self.node_value[nid] = value;
        }
    }

    pub fn transfer_quad_node(&mut self, dom: &Domain, itg: &IntegralDomain) {
        // skip if not function type
        if self.scl_type != ScalarDomainType::Function {
            return;
        }

        // initialize global matrix and vector
        let num_node = dom.num_node;
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let mut b_vec = Col::zeros(num_node);

        // iterate over elements
        for eid in 0..dom.num_elem {
            // initialize local matrices
            let num_node_loc = dom.elem_node[eid];
            let mut m_loc = vec![vec![0.0; num_node_loc]; num_node_loc];
            let mut b_loc = vec![0.0; num_node_loc];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];

            // get properties
            let quad_val = &self.quad_value[eid];

            // assemble local mass matrix and rhs vector
            match num_node_loc {
                3 => {
                    for qid in 0..num_quad {
                        let n = tri3_eval(A_TRI3[qid], B_TRI3[qid]);
                        let coeff = W_TRI3[qid] * jac_det[qid];
                        for v in 0..num_node_loc {
                            for j in 0..num_node_loc {
                                m_loc[v][j] += coeff * n[v] * n[j];
                            }
                            b_loc[v] += coeff * quad_val[qid] * n[v];
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        let n = quad4_eval(A_QUAD4[qid], B_QUAD4[qid]);
                        let coeff = W_QUAD4[qid] * jac_det[qid];
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
            let node_id = &dom.elem_node_id[eid];
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
            .map_err(|_| FEChemError::FailedMatrixSolve {caller: "ScalarDomain::transfer_quad_node".to_string()}).unwrap();
        let symbolic = SymbolicLu::try_new(m_mat.symbolic())
            .map_err(|_| FEChemError::FailedMatrixSolve {caller: "ScalarDomain::transfer_quad_node".to_string()}).unwrap();
        let lu = Lu::try_new_with_symbolic(symbolic, m_mat.as_ref())
            .map_err(|_| FEChemError::FailedMatrixSolve {caller: "ScalarDomain::transfer_quad_node".to_string()}).unwrap();
        let u_vec = lu.solve(&b_vec);

        // store nodal values
        for nid in 0..num_node {
            self.node_value[nid] = u_vec[nid];
        }

    }

}

// this was separated from the struct to prevent borrowing issues
pub fn scldom_update_function(vars: &mut Variables, scldom_id: usize, t: f64) {
    if vars.scl_dom[scldom_id].scl_type != ScalarDomainType::Function {
        return;
    }

    // get objects
    let dom_id = vars.scl_dom[scldom_id].dom_id;
    let dom = &vars.dom[dom_id];
    let itgdom = &vars.itg_dom[dom_id];

    // iterate over quadrature points
    for eid in 0..dom.num_elem {
        for qid in 0..itgdom.num_quad[eid] {
            let value = scldom_quad_value(vars, scldom_id, eid, qid, t);
            vars.scl_dom[scldom_id].quad_value[eid][qid] = value;
        }
    }
}

pub fn scldom_quad_value(vars: &Variables, scldom_id: usize, eid: usize, qid: usize, t: f64) -> f64 {
    let scldom = &vars.scl_dom[scldom_id];
    let dom = &vars.dom[scldom.dom_id];
    let itgdom = &vars.itg_dom[scldom.dom_id];
    let x = itgdom.quad_x[eid][qid];
    let y = itgdom.quad_y[eid][qid];

    // dependencies are unknown-type scalars only
    let num_node_loc = dom.elem_node[eid];
    let node_id = &dom.elem_node_id[eid];
    let mut scl_val: Vec<f64> = Vec::new();
    for &scldom_sub in &scldom.fun_sclid {
        let node_value = &vars.scl_dom[scldom_sub].node_value;
        let val = match num_node_loc {
            3 => {
                let n = tri3_eval(A_TRI3[qid], B_TRI3[qid]);
                let mut val = 0.0;
                for v in 0..num_node_loc {
                    val += n[v] * node_value[node_id[v]];
                }
                val
            }
            4 => {
                let n = quad4_eval(A_QUAD4[qid], B_QUAD4[qid]);
                let mut val = 0.0;
                for v in 0..num_node_loc {
                    val += n[v] * node_value[node_id[v]];
                }
                val
            }
            _ => { panic!("Invalid element type"); }
        };
        scl_val.push(val);
    }
    (scldom.fun_func)(t, [x, y], scl_val.as_slice())
}
