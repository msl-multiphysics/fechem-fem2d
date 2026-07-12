use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::vars::Variables;
use crate::base::write_csv::write_vecbnd_csv;
use crate::base::write_vtu::write_vecbnd_vtu;

pub enum VectorBoundaryType {
    Constant {
        value_x: f64,
        value_y: f64,
    },
    Function {
        func: Box<dyn Fn(f64, &[f64]) -> (f64, f64) + Send + Sync>,    // function of unknown scalars only
        scldom_ids: Vec<usize>,
    },
}

impl Default for VectorBoundaryType {
    fn default() -> VectorBoundaryType {
        VectorBoundaryType::Constant { value_x: 0.0, value_y: 0.0 }
    }
}

#[derive(Default)]
pub struct VectorBoundary {
    // ids
    pub vecbnd_id: usize,
    pub bnd_id: usize, // boundary this scalar is attached to

    // values
    pub vec_type: VectorBoundaryType,
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied
    pub node_value_x: Vec<f64>, // [nid] -> x values at nodes
    pub node_value_y: Vec<f64>, // [nid] -> y values at nodes

    // output file
    pub file_name: String, // path to file without extension
    pub file_type: String, // leave empty if no output
}

impl VectorBoundary {
    pub fn new_from_constant(vecbnd_id: usize, bnd: &Boundary, value_const_x: f64, value_const_y: f64, file_path: String) -> Result<VectorBoundary, FEChemError> {
        // create struct
        let mut vecbnd = VectorBoundary::default();
        vecbnd.vecbnd_id = vecbnd_id;
        vecbnd.bnd_id = bnd.bnd_id;

        // set values
        vecbnd.vec_type = VectorBoundaryType::Constant { value_x: value_const_x, value_y: value_const_y };
        vecbnd.node_dir = vec![false; bnd.num_node];
        vecbnd.node_value_x = vec![value_const_x; bnd.num_node];
        vecbnd.node_value_y = vec![value_const_y; bnd.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            vecbnd.file_name = String::new();
            vecbnd.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            if parts.len() < 2 {
                return Err(FEChemError::InvalidOutputPath {
                    caller: "VectorBoundary::new_from_constant".to_string(),
                    file_path,
                });
            }
            vecbnd.file_name = parts[0..parts.len() - 1].join(".");
            vecbnd.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(vecbnd)
    }

    pub fn new_from_function(vecbnd_id: usize, bnd: &Boundary, value_func: Box<dyn Fn(f64, &[f64]) -> (f64, f64) + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<VectorBoundary, FEChemError> {
        // create struct
        let mut vecbnd = VectorBoundary::default();
        vecbnd.vecbnd_id = vecbnd_id;
        vecbnd.bnd_id = bnd.bnd_id;

        // non-constant properties are evaluated at quadrature points
        // these will be updated only when writing
        vecbnd.vec_type = VectorBoundaryType::Function {
            func: value_func,
            scldom_ids: scldom_ids,
        };
        vecbnd.node_dir = vec![false; bnd.num_node];
        vecbnd.node_value_x = vec![0.0; bnd.num_node];
        vecbnd.node_value_y = vec![0.0; bnd.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            vecbnd.file_name = String::new();
            vecbnd.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            if parts.len() < 2 {
                return Err(FEChemError::InvalidOutputPath {
                    caller: "VectorBoundary::new_from_function".to_string(),
                    file_path,
                });
            }
            vecbnd.file_name = parts[0..parts.len() - 1].join(".");
            vecbnd.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(vecbnd)
    }

    pub fn write(&self, bnd: &Boundary, ts: usize) -> Result<(), FEChemError> {
        // // write depending on file type
        match self.file_type.as_str() {
            "csv" => write_vecbnd_csv(&bnd, &self, ts)?,
            "vtu" => write_vecbnd_vtu(&bnd, &self, ts)?,
            "" => (), // do nothing if file type is empty
            _ => {
                return Err(FEChemError::UnsupportedFileFormat {
                    caller: "VectorBoundary::write".to_string(),
                    type_need: "csv or vtu".to_string(),
                    type_got: self.file_type.clone(),
                });
            }
        }

        // placeholder
        Ok(())
    }

    pub fn compute_quad(&self, vars: &Variables, eid: usize, qid: usize, t: f64) -> (f64, f64) {
        // evaluate based on type
        match &self.vec_type {
            VectorBoundaryType::Constant { value_x, value_y } => {
                return (*value_x, *value_y);
            }
            VectorBoundaryType::Function { func, scldom_ids } => {
                // get boundary
                let bnd = &vars.bnd[self.bnd_id];
                let itgbnd = &vars.itg_bnd[self.bnd_id];

                // get scalar values
                let mut val = Vec::new();
                for &scldom_id in scldom_ids {
                    let scldom_sub = &vars.scl_dom[scldom_id];
                    let val_sub = scldom_sub.compute_quad_unknown_boundary(bnd, itgbnd, eid, qid);
                    val.push(val_sub);
                }

                // evaluate function
                return func(t, &val);
            }
        }
    }
}

// put outside of struct to avoid circular dependency
pub fn update_function_vecbnd(vars: &mut Variables, vecbnd_id: usize, t: f64) {
    // skip if not function type
    let vecbnd = &vars.vec_bnd[vecbnd_id]; // immutable borrow
    let (func, scldom_ids) = match &vecbnd.vec_type {
        VectorBoundaryType::Function { func, scldom_ids } => (func, scldom_ids),
        _ => {
            return;
        }
    };

    // initialize value vector
    let bnd = &vars.bnd[vecbnd.bnd_id];
    let mut node_value_x = Vec::with_capacity(bnd.num_node);
    let mut node_value_y = Vec::with_capacity(bnd.num_node);

    // evaluate function on node points
    for nid in 0..bnd.num_node {
        let mut val = Vec::new();
        for &scldom_id in scldom_ids {
            let scldom_sub = &vars.scl_dom[scldom_id];
            let nid_dom = bnd.node_bnd_dom_id[nid];
            let val_sub = scldom_sub.node_value[nid_dom];
            val.push(val_sub);
        }
        let (value_x, value_y) = func(t, &val);
        node_value_x.push(value_x);
        node_value_y.push(value_y);
    }

    // update node values
    let vecbnd = &mut vars.vec_bnd[vecbnd_id]; // mutable borrow
    vecbnd.node_value_x = node_value_x;
    vecbnd.node_value_y = node_value_y;
}
