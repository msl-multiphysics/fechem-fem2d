use crate::base::error::FEChemError;
use crate::base::geom_itf::Interface;
use faer::Col;

pub enum VectorInterfaceType {
    Constant { value_x: f64, value_y: f64 },
    Unknown { start: usize },
}

impl Default for VectorInterfaceType {
    fn default() -> VectorInterfaceType {
        VectorInterfaceType::Constant { value_x: 0.0, value_y: 0.0 }
    }
}

#[derive(Default)]
pub struct VectorInterface {
    // ids
    pub vecitf_id: usize,
    pub itf_id: usize, // interface this vector is attached to

    // values
    pub vec_type: VectorInterfaceType,
    pub node_value_x: Vec<f64>, // [nid] -> x values at nodes
    pub node_value_y: Vec<f64>, // [nid] -> y values at nodes
    pub node_prev_x: Vec<f64>,  // [nid] -> x values at previous time step
    pub node_prev_y: Vec<f64>,  // [nid] -> y values at previous time step
    pub node_dir: Vec<bool>,  // [nid] -> true if dirichlet BC is applied

    // output file
    pub file_name: String, // path to file without extension
    pub file_type: String, // leave empty if no output
}

impl VectorInterface {
    pub fn new_from_constant(vecitf_id: usize, itf: &Interface, value_const_x: f64, value_const_y: f64, file_path: String) -> Result<VectorInterface, FEChemError> {
        // create struct
        let mut vecitf = VectorInterface::default();
        vecitf.vecitf_id = vecitf_id;
        vecitf.itf_id = itf.itf_id;

        // set values
        vecitf.vec_type = VectorInterfaceType::Constant { value_x: value_const_x, value_y: value_const_y };
        vecitf.node_value_x = vec![value_const_x; itf.num_node];
        vecitf.node_value_y = vec![value_const_y; itf.num_node];
        vecitf.node_prev_x = vec![value_const_x; itf.num_node];
        vecitf.node_prev_y = vec![value_const_y; itf.num_node];
        vecitf.node_dir = vec![false; itf.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            vecitf.file_name = String::new();
            vecitf.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            vecitf.file_name = parts[0..parts.len() - 1].join(".");
            vecitf.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(vecitf)
    }

    // TODO: implement new_from_function

    pub fn new_from_unknown(vecitf_id: usize, itf: &Interface, value_init_x: f64, value_init_y: f64, file_path: String) -> Result<VectorInterface, FEChemError> {
        // create struct
        let mut vecitf = VectorInterface::default();
        vecitf.vecitf_id = vecitf_id;
        vecitf.itf_id = itf.itf_id;

        // set values
        vecitf.vec_type = VectorInterfaceType::Unknown { start: 0 };
        vecitf.node_value_x = vec![value_init_x; itf.num_node];
        vecitf.node_value_y = vec![value_init_y; itf.num_node];
        vecitf.node_prev_x = vec![value_init_x; itf.num_node];
        vecitf.node_prev_y = vec![value_init_y; itf.num_node];
        vecitf.node_dir = vec![false; itf.num_node];

        // set outputs if file path is not empty
        if file_path == "" {
            vecitf.file_name = String::new();
            vecitf.file_type = String::new();
        } else {
            let parts: Vec<&str> = file_path.split('.').collect();
            vecitf.file_name = parts[0..parts.len() - 1].join(".");
            vecitf.file_type = parts[parts.len() - 1].to_string();
        }

        // result
        Ok(vecitf)
    }

    // TODO: implement write

    pub fn update_unknown(&mut self, itf: &Interface, x_vec: &Col<f64>) {
        match self.vec_type {
            VectorInterfaceType::Unknown { start } => {
                for nid in 0..itf.num_node {
                    let xid_x = start + nid;
                    let xid_y = start + nid + itf.num_node;
                    self.node_value_x[nid] = x_vec[xid_x];
                    self.node_value_y[nid] = x_vec[xid_y];
                }
            }
            _ => {
                return;
            }
        }
    }

    pub fn update_prev(&mut self) {
        match self.vec_type {
            VectorInterfaceType::Unknown { .. } => {
                self.node_prev_x.copy_from_slice(&self.node_value_x);
                self.node_prev_y.copy_from_slice(&self.node_value_y);
            }
            _ => {
                return;
            }
        }
    }
}
