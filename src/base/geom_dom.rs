use crate::base::error::FEChemError;
use crate::base::mesh::Mesh;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct Domain {
    // ids
    pub dom_id: usize,

    // node data
    pub num_node: usize,
    pub node_x: Vec<f64>,
    pub node_y: Vec<f64>,

    // element data
    pub num_elem: usize,  // number of 2d elements
    pub elem_node_num: Vec<usize>,  // number of nodes (3 - tri; 4 - quad)
    pub elem_node_id: Vec<Vec<usize>>,  // node ids of elements

    // domain-mesh mapping
    pub node_dom_mesh_id: Vec<usize>,
    pub node_mesh_dom_id: HashMap<usize, usize>,
    pub elem_dom_mesh_id: Vec<usize>,
    pub elem_mesh_dom_id: HashMap<usize, usize>,

}

impl Domain {
    pub fn new(dom_id: usize, mesh: &Mesh, reg_id: usize) -> Result<Domain, FEChemError> {
        // initialize domain
        let mut dom = Domain::default();
        dom.dom_id = dom_id;

        // step 1: build domain-mesh mappings

        // create domain-mesh element mapping
        dom.elem_dom_mesh_id = mesh.reg2d_elem_id[reg_id].clone();
        for (elem_did, &elem_mid) in dom.elem_dom_mesh_id.iter().enumerate() {
            dom.elem_mesh_dom_id.insert(elem_mid, elem_did);
        }

        // build list of mesh node ids
        let mut node_set = HashSet::new();  
        for &elem_mid in dom.elem_dom_mesh_id.iter() {
            for &node_mid in &mesh.elm2d_node_id[elem_mid] {
                node_set.insert(node_mid);
            }
        }
        
        // create domain-mesh node mapping
        dom.node_dom_mesh_id = node_set.into_iter().collect::<Vec<usize>>();
        for (node_did, &node_mid) in dom.node_dom_mesh_id.iter().enumerate() {
            dom.node_mesh_dom_id.insert(node_mid, node_did);
        }

        // step 2: get data from mesh

        // get node coordinates
        for &node_mid in dom.node_dom_mesh_id.iter() {
            dom.node_x.push(mesh.node_x[node_mid]);
            dom.node_y.push(mesh.node_y[node_mid]);
        }
        dom.num_node = dom.node_dom_mesh_id.len();

        // get element data
        for &elem_mid in dom.elem_dom_mesh_id.iter() {
            // element type - copy from mesh
            dom.elem_node_num.push(mesh.elm2d_node_num[elem_mid]);

            // node ids - convert from global to local
            let mut node_id = Vec::with_capacity(mesh.elm2d_node_num[elem_mid]);
            for (_, &node_mid) in mesh.elm2d_node_id[elem_mid].iter().enumerate() {
                node_id.push(dom.node_mesh_dom_id[&node_mid]);
            }
            dom.elem_node_id.push(node_id);
        }
        dom.num_elem = dom.elem_dom_mesh_id.len();

        // result
        Ok(dom)

    }
}
