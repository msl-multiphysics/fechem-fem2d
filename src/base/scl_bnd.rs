use crate::base::geom_bnd::Boundary;
use crate::base::error::FEChemError;
use crate::base::write_csv::write_sclbnd_csv;
use crate::base::write_vtu::write_sclbnd_vtu;

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
}
