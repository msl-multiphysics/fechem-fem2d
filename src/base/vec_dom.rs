use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::itg_bnd::IntegralBoundary;
use crate::base::itg_dom::IntegralDomain;
use crate::base::vars::Variables;
use crate::base::write_csv::write_vecdom_csv;
use crate::base::write_vtu::write_vecdom_vtu;
use faer::Col;

pub enum VectorDomainType {
    Constant {
        value_x: f64,
        value_y: f64,
    },
    Function {
        func: Box<dyn Fn(f64, &[f64]) -> (f64, f64) + Send + Sync>,  // function of unknown scalars only
        scldom_ids: Vec<usize>,
    },
    Unknown {
        start: usize,
    },
}

impl Default for VectorDomainType {
    fn default() -> VectorDomainType {
        VectorDomainType::Constant { value_x: 0.0, value_y: 0.0 }
    }
}

#[derive(Default)]
pub struct VectorDomain {
    // ids
    pub vecdom_id: usize,
    pub dom_id: usize, // domain this scalar is attached to

    // values
    pub vec_type: VectorDomainType,
    pub node_dir: Vec<bool>,          // [nid] -> true if dirichlet BC is applied
    pub node_value_x: Vec<f64>,       // [nid] -> x values at nodes
    pub node_value_y: Vec<f64>,       // [nid] -> y values at nodes
    pub node_value_prev_x: Vec<f64>,  // [nid] -> x values at previous time step
    pub node_value_prev_y: Vec<f64>,  // [nid] -> y values at previous time step

    // output file
    pub file_name: String, // path to file without extension
    pub file_type: String, // leave empty if no output
}

impl VectorDomain {
    pub fn new_from_constant(vecdom_id: usize, dom: &Domain, value_const_x: f64, value_const_y: f64, file_path: String) -> Result<VectorDomain, FEChemError> {
        // create struct
        let mut vecdom = VectorDomain::default();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // set values
        vecdom.vec_type = VectorDomainType::Constant { value_x: value_const_x, value_y: value_const_y };
        vecdom.node_dir = vec![false; dom.num_node];
        vecdom.node_value_x = vec![value_const_x; dom.num_node];
        vecdom.node_value_y = vec![value_const_y; dom.num_node];
        vecdom.node_value_prev_x = vec![value_const_x; dom.num_node];
        vecdom.node_value_prev_y = vec![value_const_y; dom.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            vecdom.file_name = String::new();
            vecdom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            vecdom.file_name = parts[0..parts.len() - 1].join(".");
            vecdom.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(vecdom)
    }

    pub fn new_from_function(vecdom_id: usize, dom: &Domain, value_func: Box<dyn Fn(f64, &[f64]) -> (f64, f64) + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<VectorDomain, FEChemError> {
        // create struct
        let mut vecdom = VectorDomain::default();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // non-constant properties are evaluated at quadrature points
        // these will be updated only when writing
        vecdom.vec_type = VectorDomainType::Function {
            func: value_func,
            scldom_ids: scldom_ids,
        };
        vecdom.node_dir = vec![false; dom.num_node];
        vecdom.node_value_x = vec![0.0; dom.num_node];
        vecdom.node_value_y = vec![0.0; dom.num_node];
        vecdom.node_value_prev_x = vec![0.0; dom.num_node];
        vecdom.node_value_prev_y = vec![0.0; dom.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            vecdom.file_name = String::new();
            vecdom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            vecdom.file_name = parts[0..parts.len() - 1].join(".");
            vecdom.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(vecdom)
    }

    pub fn new_from_unknown(vecdom_id: usize, dom: &Domain, value_init_x: f64, value_init_y: f64, file_path: String) -> Result<VectorDomain, FEChemError> {
        // create struct
        let mut vecdom = VectorDomain::default();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // set values
        vecdom.vec_type = VectorDomainType::Unknown { start: 0 };
        vecdom.node_dir = vec![false; dom.num_node];
        vecdom.node_value_x = vec![value_init_x; dom.num_node];
        vecdom.node_value_y = vec![value_init_y; dom.num_node];
        vecdom.node_value_prev_x = vec![value_init_x; dom.num_node];
        vecdom.node_value_prev_y = vec![value_init_y; dom.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            vecdom.file_name = String::new();
            vecdom.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            vecdom.file_name = parts[0..parts.len() - 1].join(".");
            vecdom.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(vecdom)
    }

    pub fn write(&self, dom: &Domain, ts: usize) -> Result<(), FEChemError> {
        // TODO: update if function type

        // write depending on file type
        match self.file_type.as_str() {
            "csv" => write_vecdom_csv(&dom, &self, ts)?,
            "vtu" => write_vecdom_vtu(&dom, &self, ts)?,
            "" => (), // do nothing if file type is empty
            _ => panic!("Unsupported file type: {}", self.file_type),
        }

        // placeholder
        Ok(())
    }

    pub fn compute_quad(&self, vars: &Variables, eid: usize, qid: usize, t: f64) -> (f64, f64) {
        // evaluate based on type
        match &self.vec_type {
            VectorDomainType::Constant { value_x, value_y } => {
                return (*value_x, *value_y);
            }
            VectorDomainType::Unknown { .. } => {
                let dom = &vars.dom[self.dom_id];
                let itgdom = &vars.itg_dom[self.dom_id];
                return self.compute_quad_unknown_domain(dom, itgdom, eid, qid);
            }
            VectorDomainType::Function { func, scldom_ids } => {
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

    pub fn compute_quad_unknown_domain(&self, dom: &Domain, itgdom: &IntegralDomain, eid: usize, qid: usize) -> (f64, f64) {
        let num_node = dom.elem_node[eid];
        let mut val_x = 0.0;
        let mut val_y = 0.0;
        for v in 0..num_node {
            let nid = dom.elem_node_id[eid][v];
            val_x += itgdom.quad_n[eid][qid][v] * self.node_value_x[nid];
            val_y += itgdom.quad_n[eid][qid][v] * self.node_value_y[nid];
        }
        return (val_x, val_y);
    }

    pub fn compute_quad_grad_unknown_domain(&self, dom: &Domain, itgdom: &IntegralDomain, eid: usize, qid: usize) -> [[f64; 2]; 2] {
        let num_node = dom.elem_node[eid];
        let mut grad_x = [0.0, 0.0];
        let mut grad_y = [0.0, 0.0];
        for v in 0..num_node {
            let nid = dom.elem_node_id[eid][v];
            let gnx = itgdom.quad_gnx[eid][qid][v];
            let gny = itgdom.quad_gny[eid][qid][v];
            grad_x[0] += gnx * self.node_value_x[nid];
            grad_x[1] += gny * self.node_value_x[nid];
            grad_y[0] += gnx * self.node_value_y[nid];
            grad_y[1] += gny * self.node_value_y[nid];
        }
        [grad_x, grad_y]
    }

    pub fn compute_quad_unknown_boundary(&self, bnd: &Boundary, itgbnd: &IntegralBoundary, eid: usize, qid: usize) -> (f64, f64) {
        let num_node = bnd.elem_node[eid];
        let mut val_x = 0.0;
        let mut val_y = 0.0;
        for v in 0..num_node {
            let nid_bnd = bnd.elem_node_id[eid][v];
            let nid_dom = bnd.node_bnd_dom_id[nid_bnd];
            val_x += itgbnd.quad_n[eid][qid][v] * self.node_value_x[nid_dom];
            val_y += itgbnd.quad_n[eid][qid][v] * self.node_value_y[nid_dom];
        }
        return (val_x, val_y);
    }

    pub fn compute_quad_prev(&self, vars: &Variables, eid: usize, qid: usize, t_prev: f64) -> (f64, f64) {
        // evaluate based on type
        match &self.vec_type {
            VectorDomainType::Constant { value_x, value_y } => {
                return (*value_x, *value_y);
            }
            VectorDomainType::Unknown { .. } => {
                let dom = &vars.dom[self.dom_id];
                let itgdom = &vars.itg_dom[self.dom_id];
                return self.compute_quad_unknown_domain_prev(dom, itgdom, eid, qid);
            }
            VectorDomainType::Function { func, scldom_ids } => {
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

    pub fn compute_quad_unknown_domain_prev(&self, dom: &Domain, itgdom: &IntegralDomain, eid: usize, qid: usize) -> (f64, f64) {
        let num_node = dom.elem_node[eid];
        let mut val_x = 0.0;
        let mut val_y = 0.0;
        for v in 0..num_node {
            let nid = dom.elem_node_id[eid][v];
            val_x += itgdom.quad_n[eid][qid][v] * self.node_value_prev_x[nid];
            val_y += itgdom.quad_n[eid][qid][v] * self.node_value_prev_y[nid];
        }
        return (val_x, val_y);
    }

    pub fn compute_quad_grad_unknown_domain_prev(&self, dom: &Domain, itgdom: &IntegralDomain, eid: usize, qid: usize) -> [[f64; 2]; 2] {
        let num_node = dom.elem_node[eid];
        let mut grad_x = [0.0, 0.0];
        let mut grad_y = [0.0, 0.0];
        for v in 0..num_node {
            let nid = dom.elem_node_id[eid][v];
            let gnx = itgdom.quad_gnx[eid][qid][v];
            let gny = itgdom.quad_gny[eid][qid][v];
            grad_x[0] += gnx * self.node_value_prev_x[nid];
            grad_x[1] += gny * self.node_value_prev_x[nid];
            grad_y[0] += gnx * self.node_value_prev_y[nid];
            grad_y[1] += gny * self.node_value_prev_y[nid];
        }
        [grad_x, grad_y]
    }

    pub fn update_unknown(&mut self, dom: &Domain, x_vec: &Col<f64>) {
        // skip if not unknown type
        match self.vec_type {
            VectorDomainType::Unknown { start } => {
                let num_node = dom.num_node;
                for nid in 0..num_node {
                    let xid_x = start + nid;
                    let xid_y = start + nid + dom.num_node;
                    let value_x = x_vec[xid_x];
                    let value_y = x_vec[xid_y];
                    self.node_value_x[nid] = value_x;
                    self.node_value_y[nid] = value_y;
                }
            }
            _ => {
                return;
            }
        }
    }

    pub fn update_prev(&mut self) {
        match self.vec_type {
            VectorDomainType::Unknown { .. } => {
                self.node_value_prev_x.copy_from_slice(&self.node_value_x);
                self.node_value_prev_y.copy_from_slice(&self.node_value_y);
            }
            _ => {
                return;
            }
        }
    }

}

// put outside of struct to avoid circular dependency
pub fn update_function_vecdom(vars: &mut Variables, vecdom_id: usize, t: f64) {
    // skip if not function type
    let vecdom = &vars.vec_dom[vecdom_id]; // immutable borrow
    let (func, scldom_ids) = match &vecdom.vec_type {
        VectorDomainType::Function { func, scldom_ids } => (func, scldom_ids),
        _ => { return; }
    };

    // initialize value vector
    let dom = &vars.dom[vecdom.dom_id];
    let mut node_value_x = Vec::with_capacity(dom.num_node);
    let mut node_value_y = Vec::with_capacity(dom.num_node);

    // evaluate function on node points
    for nid in 0..dom.num_node {
        let mut val = Vec::new();
        for &scldom_id in scldom_ids {
            let scldom_sub = &vars.scl_dom[scldom_id];
            let val_sub = scldom_sub.node_value[nid];
            val.push(val_sub);
        }
        let (value_x, value_y) = func(t, &val);
        node_value_x.push(value_x);
        node_value_y.push(value_y);
    }

    // update node values
    let vecdom = &mut vars.vec_dom[vecdom_id]; // mutable borrow
    vecdom.node_value_x = node_value_x;
    vecdom.node_value_y = node_value_y;
}
