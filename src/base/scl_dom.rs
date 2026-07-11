use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::itg_bnd::IntegralBoundary;
use crate::base::itg_dom::IntegralDomain;
use crate::base::vars::Variables;
use crate::base::write_csv::write_scldom_csv;
use crate::base::write_vtu::write_scldom_vtu;
use faer::Col;

pub enum ScalarDomainType {
    Constant {
        value: f64,
    },
    Function {
        func: Box<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>,  // function of unknown scalars only
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
    pub node_dir: Vec<bool>,        // [nid] -> true if dirichlet BC is applied
    pub node_value: Vec<f64>,       // [nid] -> values at nodes
    pub node_value_prev: Vec<f64>,  // [nid] -> values at previous time step

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
        scldom.node_dir = vec![false; dom.num_node];
        scldom.node_value = vec![value_const; dom.num_node];
        scldom.node_value_prev = vec![value_const; dom.num_node];
        
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

    pub fn new_from_function(scldom_id: usize, dom: &Domain, value_func: Box<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<ScalarDomain, FEChemError> {
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
        scldom.node_dir = vec![false; dom.num_node];
        scldom.node_value = vec![0.0; dom.num_node];
        scldom.node_value_prev = vec![0.0; dom.num_node];

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
        scldom.node_dir = vec![false; dom.num_node];
        scldom.node_value = vec![value_init; dom.num_node];
        scldom.node_value_prev = vec![value_init; dom.num_node];

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
                let itgdom = &vars.itg_dom[self.dom_id];
                return self.compute_quad_unknown_domain(dom, itgdom, eid, qid);
            }
            ScalarDomainType::Function { func, scldom_ids } => {
                // get domain
                let dom = &vars.dom[self.dom_id];
                let itgdom = &vars.itg_dom[self.dom_id];

                // get scalar values
                let mut val = Vec::new();
                for &scldom_id in scldom_ids {
                    let scldom_sub = &vars.scl_dom[scldom_id];
                    let val_sub = scldom_sub.compute_quad_unknown_domain(dom, itgdom, eid, qid);
                    val.push(val_sub);
                }

                // evaluate function
                return func(t, &val);
            }
        }
    }

    pub fn compute_quad_unknown_domain(&self, dom: &Domain, itgdom: &IntegralDomain, eid: usize, qid: usize) -> f64 {
        let num_node = dom.elem_node[eid];
        let mut value_quad = 0.0;
        for v in 0..num_node {
            let nid = dom.elem_node_id[eid][v];
            let value_node = self.node_value[nid];
            value_quad += itgdom.quad_n[eid][qid][v] * value_node;
        }
        return value_quad;
    }

    pub fn compute_quad_unknown_boundary(&self, bnd: &Boundary, itgbnd: &IntegralBoundary, eid: usize, qid: usize) -> f64 {
        let num_node = bnd.elem_node[eid];
        let mut value_quad = 0.0;
        for v in 0..num_node {
            let nid_bnd = bnd.elem_node_id[eid][v];
            let nid_dom = bnd.node_bnd_dom_id[nid_bnd];
            let value_node = self.node_value[nid_dom];
            value_quad += itgbnd.quad_n[eid][qid][v] * value_node;
        }
        return value_quad;
    }

    pub fn compute_quad_prev(&self, vars: &Variables, eid: usize, qid: usize, t_prev: f64) -> f64 {
        // evaluate based on type
        match &self.scl_type {
            ScalarDomainType::Constant { value } => {
                return *value;
            }
            ScalarDomainType::Unknown { .. } => {
                let dom = &vars.dom[self.dom_id];
                let itgdom = &vars.itg_dom[self.dom_id];
                return self.compute_quad_unknown_domain_prev(dom, itgdom, eid, qid);
            }
            ScalarDomainType::Function { func, scldom_ids } => {
                // get domain
                let dom = &vars.dom[self.dom_id];
                let itgdom = &vars.itg_dom[self.dom_id];

                // get scalar values
                let mut val = Vec::new();
                for &scldom_id in scldom_ids {
                    let scldom_sub = &vars.scl_dom[scldom_id];
                    let val_sub = scldom_sub.compute_quad_unknown_domain_prev(dom, itgdom, eid, qid);
                    val.push(val_sub);
                }

                // evaluate function
                return func(t_prev, &val);
            }
        }
    }

    pub fn compute_quad_unknown_domain_prev(&self, dom: &Domain, itgdom: &IntegralDomain, eid: usize, qid: usize) -> f64 {
        let num_node = dom.elem_node[eid];
        let mut value_quad = 0.0;
        for v in 0..num_node {
            let nid = dom.elem_node_id[eid][v];
            let value_node = self.node_value_prev[nid];
            value_quad += itgdom.quad_n[eid][qid][v] * value_node;
        }
        return value_quad;
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
                self.node_value_prev.copy_from_slice(&self.node_value);
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
        let mut val = Vec::new();
        for &scldom_id in scldom_ids {
            let scldom_sub = &vars.scl_dom[scldom_id];
            let val_sub = scldom_sub.node_value[nid];
            val.push(val_sub);
        }
        let value = func(t, &val);
        node_value.push(value);
    }

    // update node values
    let scldom = &mut vars.scl_dom[scldom_id]; // mutable borrow
    scldom.node_value = node_value;
}
