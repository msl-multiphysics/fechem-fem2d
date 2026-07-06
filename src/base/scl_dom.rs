use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::vars::Variables;
use crate::base::write_csv::write_scldom_csv;
use crate::base::write_vtu::write_scldom_vtu;
use crate::shape::prelude::*;
use faer::Col;

pub enum ScalarDomainType {
    Constant {
        value: f64,
    },
    Function {
        func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64 + Send + Sync>,  // function of unknown scalars only
        scldom_ids: Vec<usize>,
    },
    Unknown {
        start: usize,
    },
}

impl Default for ScalarDomainType {
    fn default() -> ScalarDomainType {
        ScalarDomainType::Constant { value: 0.0 }
    }
}

#[derive(Default)]
pub struct ScalarDomain {
    // ids
    pub scldom_id: usize,
    pub dom_id: usize, // domain this scalar is attached to

    // values
    pub scl_type: ScalarDomainType,
    pub node_value: Vec<f64>, // [nid] -> values at nodes
    pub node_prev: Vec<f64>,  // [nid] -> values at previous time step
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied

    // output file
    pub file_name: String, // path to file without extension
    pub file_type: String, // leave empty if no output
}

impl ScalarDomain {
    pub fn new_from_constant(scldom_id: usize, dom: &Domain, value_const: f64, file_path: String) -> Result<ScalarDomain, FEChemError> {
        // create struct
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // set values
        scldom.scl_type = ScalarDomainType::Constant { value: value_const };
        scldom.node_value = vec![value_const; dom.num_node];
        scldom.node_prev = vec![value_const; dom.num_node];
        scldom.node_dir = vec![false; dom.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            scldom.file_name = parts[0..parts.len() - 1].join(".");
            scldom.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(scldom)
    }

    pub fn new_from_function(scldom_id: usize, dom: &Domain, value_func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64 + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<ScalarDomain, FEChemError> {
        // create struct
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // non-constant properties are evaluated at quadrature points
        // these will be updated only when writing
        scldom.scl_type = ScalarDomainType::Function {
            func: value_func,
            scldom_ids: scldom_ids,
        };
        scldom.node_value = vec![0.0; dom.num_node];
        scldom.node_prev = vec![0.0; dom.num_node];
        scldom.node_dir = vec![false; dom.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            scldom.file_name = parts[0..parts.len() - 1].join(".");
            scldom.file_type = parts[parts.len() - 1].to_string();
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
        scldom.scl_type = ScalarDomainType::Unknown { start: 0 };
        scldom.node_value = vec![value_init; dom.num_node];
        scldom.node_prev = vec![value_init; dom.num_node];
        scldom.node_dir = vec![false; dom.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            scldom.file_name = parts[0..parts.len() - 1].join(".");
            scldom.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(scldom)
    }

    pub fn write(&self, dom: &Domain, ts: usize) -> Result<(), FEChemError> {
        // TODO: update if function type

        // write depending on file type
        match self.file_type.as_str() {
            "csv" => write_scldom_csv(&dom, &self, ts)?,
            "vtu" => write_scldom_vtu(&dom, &self, ts)?,
            "" => (), // do nothing if file type is empty
            _ => panic!("Unsupported file type: {}", self.file_type),
        }

        // placeholder
        Ok(())
    }

    pub fn compute_quad(&self, vars: &Variables, eid: usize, qid: usize, t: f64) -> f64 {
        // evaluate based on type
        match &self.scl_type {
            ScalarDomainType::Constant { value } => {
                return *value;
            }
            ScalarDomainType::Unknown { .. } => {
                let dom = &vars.dom[self.dom_id];
                return self.compute_quad_unknown_domain(dom, eid, qid);
            }
            ScalarDomainType::Function { func, scldom_ids } => {
                // get coordinates
                let dom = &vars.dom[self.dom_id];
                let itgdom = &vars.itg_dom[self.dom_id];
                let x = itgdom.quad_x[eid][qid];
                let y = itgdom.quad_y[eid][qid];

                // get scalar values
                let mut val = Vec::new();
                for &scldom_id in scldom_ids {
                    let scldom_sub = &vars.scl_dom[scldom_id];
                    let val_sub = scldom_sub.compute_quad_unknown_domain(dom, eid, qid);
                    val.push(val_sub);
                }

                // evaluate function
                return func(t, [x, y], &val);
            }
        }
    }

    pub fn compute_quad_unknown_domain(&self, dom: &Domain, eid: usize, qid: usize) -> f64 {
        let num_node = dom.elem_node_num[eid];
        match num_node {
            3 => {
                let n = tri3_eval(A_TRI3[qid], B_TRI3[qid]);
                let mut val = 0.0;
                for v in 0..num_node {
                    let nid = dom.elem_node_id[eid][v];
                    val += n[v] * self.node_value[nid];
                }
                return val;
            }
            4 => {
                let n = quad4_eval(A_QUAD4[qid], B_QUAD4[qid]);
                let mut val = 0.0;
                for v in 0..num_node {
                    let nid = dom.elem_node_id[eid][v];
                    val += n[v] * self.node_value[nid];
                }
                return val;
            }
            _ => {
                panic!("Unsupported number of nodes: {}", num_node);
            }
        }
    }

    pub fn compute_quad_unknown_boundary(&self, bnd: &Boundary, eid: usize, qid: usize) -> f64 {
        let num_node = bnd.elem_node_num[eid];
        match num_node {
            2 => {
                let n = lin2_eval(A_LIN2[qid]);
                let mut val = 0.0;
                for v in 0..num_node {
                    let nid = bnd.elem_node_id[eid][v];
                    val += n[v] * self.node_value[nid];
                }
                return val;
            }
            _ => {
                panic!("Unsupported number of nodes: {}", num_node);
            }
        }
    }

    pub fn compute_quad_prev(&self, vars: &Variables, eid: usize, qid: usize, t_prev: f64) -> f64 {
        // evaluate based on type
        match &self.scl_type {
            ScalarDomainType::Constant { value } => {
                return *value;
            }
            ScalarDomainType::Unknown { .. } => {
                let dom = &vars.dom[self.dom_id];
                return self.compute_quad_unknown_domain_prev(dom, eid, qid);
            }
            ScalarDomainType::Function { func, scldom_ids } => {
                // get coordinates
                let dom = &vars.dom[self.dom_id];
                let itgdom = &vars.itg_dom[self.dom_id];
                let x = itgdom.quad_x[eid][qid];
                let y = itgdom.quad_y[eid][qid];

                // get scalar values
                let mut val = Vec::new();
                for &scldom_id in scldom_ids {
                    let scldom_sub = &vars.scl_dom[scldom_id];
                    let val_sub = scldom_sub.compute_quad_unknown_domain_prev(dom, eid, qid);
                    val.push(val_sub);
                }

                // evaluate function
                return func(t_prev, [x, y], &val);
            }
        }
    }

    pub fn compute_quad_unknown_domain_prev(&self, dom: &Domain, eid: usize, qid: usize) -> f64 {
        let num_node = dom.elem_node_num[eid];
        match num_node {
            3 => {
                let n = tri3_eval(A_TRI3[qid], B_TRI3[qid]);
                let mut val = 0.0;
                for v in 0..num_node {
                    let nid = dom.elem_node_id[eid][v];
                    val += n[v] * self.node_prev[nid];
                }
                return val;
            }
            4 => {
                let n = quad4_eval(A_QUAD4[qid], B_QUAD4[qid]);
                let mut val = 0.0;
                for v in 0..num_node {
                    let nid = dom.elem_node_id[eid][v];
                    val += n[v] * self.node_prev[nid];
                }
                return val;
            }
            _ => {
                panic!("Unsupported number of nodes: {}", num_node);
            }
        }
    }

    pub fn compute_quad_unknown_boundary_prev(&self, bnd: &Boundary, eid: usize, qid: usize) -> f64 {
        let num_node = bnd.elem_node_num[eid];
        match num_node {
            2 => {
                let n = lin2_eval(A_LIN2[qid]);
                let mut val = 0.0;
                for v in 0..num_node {
                    let nid = bnd.elem_node_id[eid][v];
                    val += n[v] * self.node_prev[nid];
                }
                return val;
            }
            _ => {
                panic!("Unsupported number of nodes: {}", num_node);
            }
        }
    }

    pub fn update_unknown(&mut self, dom: &Domain, x_vec: &Col<f64>) {
        // skip if not unknown type
        match self.scl_type {
            ScalarDomainType::Unknown { start } => {
                let num_node = dom.num_node;
                for nid in 0..num_node {
                    let xid = start + nid;
                    let value = x_vec[xid];
                    self.node_value[nid] = value;
                }
            }
            _ => {
                return;
            }
        }
    }

    pub fn update_prev(&mut self) {
        match self.scl_type {
            ScalarDomainType::Unknown { .. } => {
                self.node_prev.copy_from_slice(&self.node_value);
            }
            _ => {
                return;
            }
        }
    }

}

// put outside of struct to avoid circular dependency
pub fn update_function_scldom(vars: &mut Variables, scldom_id: usize, t: f64) {
    // skip if not function type
    let scldom = &vars.scl_dom[scldom_id]; // immutable borrow
    let (func, scldom_ids) = match &scldom.scl_type {
        ScalarDomainType::Function { func, scldom_ids } => (func, scldom_ids),
        _ => {
            return;
        }
    };

    // initialize value vector
    let dom = &vars.dom[scldom.dom_id];
    let mut node_value = Vec::with_capacity(dom.num_node);

    // evaluate function on node points
    for nid in 0..dom.num_node {
        let x = dom.node_x[nid];
        let y = dom.node_y[nid];
        let mut val = Vec::new();
        for &scldom_id in scldom_ids {
            let scldom_sub = &vars.scl_dom[scldom_id];
            let val_sub = scldom_sub.node_value[nid];
            val.push(val_sub);
        }
        let value = func(t, [x, y], &val);
        node_value.push(value);
    }

    // update node values
    let scldom = &mut vars.scl_dom[scldom_id]; // mutable borrow
    scldom.node_value = node_value;
}
