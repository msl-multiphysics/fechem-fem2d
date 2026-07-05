use crate::base::error::FEChemError;
use crate::base::geom_dom::Domain;
use crate::base::mesh::Mesh;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct Boundary {
    // ids
    pub bnd_id: usize,
    pub dom_id: usize, // domain this boundary is attached to

    // node data
    pub num_node: usize,
    pub node_x: Vec<f64>,
    pub node_y: Vec<f64>,

    // element data
    pub num_elem: usize,               // number of 1d elements
    pub elem_node_num: Vec<usize>,     // type of elements (2 - lin)
    pub elem_node_id: Vec<Vec<usize>>, // node ids of elements

    // boundary-mesh mapping
    pub node_bnd_mesh_id: Vec<usize>,
    pub node_mesh_bnd_id: HashMap<usize, usize>,
    pub elem_bnd_mesh_id: Vec<usize>,
    pub elem_mesh_bnd_id: HashMap<usize, usize>,

    // boundary-domain mapping
    pub node_bnd_dom_id: Vec<usize>,
    pub node_dom_bnd_id: HashMap<usize, usize>,
    pub elem_bnd_dom_id: Vec<usize>, // parent domain element for each boundary element
}

impl Boundary {
    pub fn new(bnd_id: usize, mesh: &Mesh, dom: &Domain, reg_id: usize) -> Result<Boundary, FEChemError> {
        // initialize boundary
        let mut bnd = Boundary::default();
        bnd.bnd_id = bnd_id;
        bnd.dom_id = dom.dom_id;

        // step 1: build boundary-mesh mappings

        // create boundary-mesh element mapping
        bnd.elem_bnd_mesh_id = mesh.reg1d_elem_id[reg_id].clone();
        for (elem_bid, &elem_mid) in bnd.elem_bnd_mesh_id.iter().enumerate() {
            bnd.elem_mesh_bnd_id.insert(elem_mid, elem_bid);
        }

        // build list of global node ids
        let mut node_set = HashSet::new();
        for &elem_mid in bnd.elem_bnd_mesh_id.iter() {
            for &node_mid in &mesh.elm1d_node_id[elem_mid] {
                node_set.insert(node_mid);
            }
        }

        // create global-local node mapping
        bnd.node_bnd_mesh_id = node_set.into_iter().collect::<Vec<usize>>();
        for (node_bid, &node_mid) in bnd.node_bnd_mesh_id.iter().enumerate() {
            bnd.node_mesh_bnd_id.insert(node_mid, node_bid);
        }

        // step 2: get data from mesh

        // get node coordinates
        for &node_mid in bnd.node_bnd_mesh_id.iter() {
            bnd.node_x.push(mesh.node_x[node_mid]);
            bnd.node_y.push(mesh.node_y[node_mid]);
        }
        bnd.num_node = bnd.node_bnd_mesh_id.len();

        // get element data
        for &elem_mid in bnd.elem_bnd_mesh_id.iter() {
            // element type - copy from mesh
            bnd.elem_node_num.push(mesh.elm1d_node_num[elem_mid]);

            // node ids - orient to match parent domain element edge traversal
            let mesh_nid = &mesh.elm1d_node_id[elem_mid];
            let d0 = dom.node_mesh_dom_id[&mesh_nid[0]];
            let d1 = dom.node_mesh_dom_id[&mesh_nid[1]];
            let (dom_eid, dom_n0, dom_n1) = find_domain_edge(dom, d0, d1)?;
            let m0 = dom.node_dom_mesh_id[dom_n0];
            let m1 = dom.node_dom_mesh_id[dom_n1];
            bnd.elem_node_id.push(vec![bnd.node_mesh_bnd_id[&m0], bnd.node_mesh_bnd_id[&m1]]);
            bnd.elem_bnd_dom_id.push(dom_eid);
        }
        bnd.num_elem = bnd.elem_bnd_mesh_id.len();

        // step 3: build boundary-domain mappings

        // iterate through boundary nodes
        for (node_bid, &node_mid) in bnd.node_bnd_mesh_id.iter().enumerate() {
            let node_did = dom.node_mesh_dom_id[&node_mid];
            bnd.node_bnd_dom_id.push(node_did);
            bnd.node_dom_bnd_id.insert(node_did, node_bid);
        }

        // result
        Ok(bnd)
    }
}

/// Find a domain element edge matching `d0` and `d1`.
/// Returns the parent element id and node ids in CCW edge order.
fn find_domain_edge(dom: &Domain, d0: usize, d1: usize) -> Result<(usize, usize, usize), FEChemError> {
    for (eid, elem) in dom.elem_node_id.iter().enumerate() {
        let n = elem.len();
        for i in 0..n {
            let a = elem[i];
            let b = elem[(i + 1) % n];
            if (a == d0 && b == d1) || (a == d1 && b == d0) {
                return Ok((eid, a, b));
            }
        }
    }
    Err(FEChemError::BoundaryEdgeNotFound {
        node0: d0,
        node1: d1,
    })
}
