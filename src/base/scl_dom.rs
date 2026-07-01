use crate::base::error::FEChemError;
use crate::base::geom_dom::Domain;
use crate::base::io_vtu::write_scldom_vtu;

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
        // // write depending on file type
        match self.file_type.as_str() {
            "vtu" => write_scldom_vtu(&dom, &self, ts)?,
            "" => (),  // do nothing if file type is empty
            _ => panic!("Unsupported file type: {}", self.file_type),
        }

        // placeholder
        Ok(())
    }

}
