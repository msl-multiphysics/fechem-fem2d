use crate::base::error::FEChemError;
use crate::base::geom_itf::Interface;
use faer::Col;

pub enum ScalarInterfaceType {
    Constant { value: f64 },
    Unknown { start: usize },
}

impl Default for ScalarInterfaceType {
    fn default() -> ScalarInterfaceType {
        ScalarInterfaceType::Constant { value: 0.0 }
    }
}

#[derive(Default)]
pub struct ScalarInterface {
    // ids
    pub sclitf_id: usize,
    pub itf_id: usize, // interface this scalar is attached to

    // values
    pub scl_type: ScalarInterfaceType,
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied
    pub node_value: Vec<f64>, // [nid] -> values at nodes
    pub node_value_prev: Vec<f64>,  // [nid] -> values at previous time step

    // output file
    pub file_name: String, // path to file without extension
    pub file_type: String, // leave empty if no output
}

impl ScalarInterface {
    pub fn new_from_constant(sclitf_id: usize, itf: &Interface, value_const: f64, file_path: String) -> Result<ScalarInterface, FEChemError> {
        // create struct
        let mut sclitf = ScalarInterface::default();
        sclitf.sclitf_id = sclitf_id;
        sclitf.itf_id = itf.itf_id;

        // set values
        sclitf.scl_type = ScalarInterfaceType::Constant { value: value_const };
        sclitf.node_dir = vec![false; itf.num_node];
        sclitf.node_value = vec![value_const; itf.num_node];
        sclitf.node_value_prev = vec![value_const; itf.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            sclitf.file_name = String::new();
            sclitf.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            sclitf.file_name = parts[0..parts.len() - 1].join(".");
            sclitf.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(sclitf)
    }

    // TODO: implement new_from_function

    pub fn new_from_unknown(sclitf_id: usize, itf: &Interface, value_init: f64, file_path: String) -> Result<ScalarInterface, FEChemError> {
        // create struct
        let mut sclitf = ScalarInterface::default();
        sclitf.sclitf_id = sclitf_id;
        sclitf.itf_id = itf.itf_id;

        // set values
        sclitf.scl_type = ScalarInterfaceType::Unknown { start: 0 };
        sclitf.node_dir = vec![false; itf.num_node];
        sclitf.node_value = vec![value_init; itf.num_node];
        sclitf.node_value_prev = vec![value_init; itf.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            sclitf.file_name = String::new();
            sclitf.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            sclitf.file_name = parts[0..parts.len() - 1].join(".");
            sclitf.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(sclitf)
    }

    // TODO: implement write

    pub fn compute_quad(&self, _eid: usize, _qid: usize, _t: f64) -> f64 {
        // evaluate based on type
        match &self.scl_type {
            ScalarInterfaceType::Constant { value } => {
                return *value;
            }
            ScalarInterfaceType::Unknown { .. } => {
                panic!("Expected constant scalar interface type.");
            }
        }
    }

    pub fn update_unknown(&mut self, itf: &Interface, x_vec: &Col<f64>) {
        match self.scl_type {
            ScalarInterfaceType::Unknown { start } => {
                for nid in 0..itf.num_node {
                    let xid = start + nid;
                    self.node_value[nid] = x_vec[xid];
                }
            }
            _ => {
                return;
            }
        }
    }

    pub fn update_prev(&mut self) {
        match self.scl_type {
            ScalarInterfaceType::Unknown { .. } => {
                self.node_value_prev.copy_from_slice(&self.node_value);
            }
            _ => {
                return;
            }
        }
    }
}
