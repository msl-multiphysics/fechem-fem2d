use crate::base::geom_bnd::Boundary;
use crate::base::error::FEChemError;
use crate::base::vars::Variables;
use crate::base::write_csv::write_sclbnd_csv;
use crate::base::write_vtu::write_sclbnd_vtu;

pub enum ScalarBoundaryType {
    Constant {
        value: f64,
    },
    Function {
        func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64 + Send + Sync>,
        scldom_ids: Vec<usize>,
    },
}

impl Default for ScalarBoundaryType {
    fn default() -> ScalarBoundaryType {
        ScalarBoundaryType::Constant { value: 0.0 }
    }
}

#[derive(Default)]
pub struct ScalarBoundary {
    // ids
    pub sclbnd_id: usize,
    pub bnd_id: usize,  // boundary this scalar is attached to

    // values
    pub scl_type: ScalarBoundaryType,
    pub node_value: Vec<f64>,  // [nid] -> values at nodes
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied

    // output file
    pub file_name: String,  // path to file without extension
    pub file_type: String,  // leave empty if no output
}

impl ScalarBoundary {
    pub fn new_from_constant(sclbnd_id: usize, bnd: &Boundary, value_const: f64, file_path: String) -> Result<ScalarBoundary, FEChemError> {
        // create struct
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // set values
        sclbnd.scl_type = ScalarBoundaryType::Constant { value: value_const };
        sclbnd.node_value = vec![value_const; bnd.num_node];
        sclbnd.node_dir = vec![false; bnd.num_node];

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

    pub fn new_from_function(sclbnd_id: usize, bnd: &Boundary, value_func: Box<dyn Fn(f64, [f64; 2], &[f64]) -> f64 + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<ScalarBoundary, FEChemError> {
        // create struct
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // non-constant properties are evaluated at quadrature points
        // these will be updated only when writing
        sclbnd.scl_type = ScalarBoundaryType::Function {
            func: value_func,
            scldom_ids: scldom_ids,
        };
        sclbnd.node_value = vec![0.0; bnd.num_node];
        sclbnd.node_dir = vec![false; bnd.num_node];

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

    pub fn compute_quad(&self, vars: &Variables, eid: usize, qid: usize, t: f64) -> f64 {
        // evaluate based on type
        match &self.scl_type {
            ScalarBoundaryType::Constant { value } => {
                return *value;
            }
            ScalarBoundaryType::Function { func, scldom_ids } => {
                // get coordinates 
                let bnd = &vars.bnd[self.bnd_id];
                let itgbnd = &vars.itg_bnd[self.bnd_id];
                let x = itgbnd.quad_x[eid][qid];
                let y = itgbnd.quad_y[eid][qid];

                // get scalar values
                let mut val = Vec::new();
                for &scldom_id in scldom_ids {
                    let scldom_sub = &vars.scl_dom[scldom_id];
                    let val_sub = scldom_sub.compute_quad_unknown_boundary(bnd, eid, qid);
                    val.push(val_sub);
                }

                // evaluate function
                return func(t, [x, y], &val);
            }
        }
    }

}

// put outside of struct to avoid circular dependency
pub fn update_function_sclbnd(vars: &mut Variables, sclbnd_id: usize, t: f64) {
    // skip if not function type
    let sclbnd = &vars.scl_bnd[sclbnd_id];  // immutable borrow
    let (func, scldom_ids) = match &sclbnd.scl_type {
        ScalarBoundaryType::Function { func, scldom_ids } => {
            (func, scldom_ids)
        }
        _ => { return; }
    };

    // initialize value vector
    let bnd = &vars.bnd[sclbnd.bnd_id];
    let mut node_value = Vec::with_capacity(bnd.num_node);

    // evaluate function on node points
    for nid in 0..bnd.num_node {
        let x = bnd.node_x[nid];
        let y = bnd.node_y[nid];
        let mut val = Vec::new();
        for &scldom_id in scldom_ids {
            let scldom_sub = &vars.scl_dom[scldom_id];
            let nid_dom = bnd.node_bnd_dom_id[nid];
            let val_sub = scldom_sub.node_value[nid_dom];
            val.push(val_sub);
        }
        let value = func(t, [x, y], &val);
        node_value.push(value);
    }

    // update node values
    let sclbnd = &mut vars.scl_bnd[sclbnd_id];  // mutable borrow
    sclbnd.node_value = node_value;

}
