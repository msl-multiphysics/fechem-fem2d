use crate::base::error::FEChemError;
use crate::base::geom_dom::Domain;
use crate::base::mesh::Mesh;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct Interface {
    // ids
    pub itf_id: usize,
    pub dom1_id: usize, // domains that interface connects
    pub dom2_id: usize, // domains that interface connects

    // node data
    pub num_node: usize,
    pub node_x: Vec<f64>,
    pub node_y: Vec<f64>,

    // element data
    // elem_node_id[eid] is arranged in CCW order for each domain
    pub num_elem: usize,                // number of 1d elements
    pub elem_node_num: Vec<usize>,      // type of elements (2 - lin)
    pub elem_node1_id: Vec<Vec<usize>>, // node ids of elements for domain 1
    pub elem_node2_id: Vec<Vec<usize>>, // node ids of elements for domain 2

    // interface-mesh mapping
    pub node_itf_mesh_id: Vec<usize>,
    pub node_mesh_itf_id: HashMap<usize, usize>,
    pub elem_itf_mesh_id: Vec<usize>,
    pub elem_mesh_itf_id: HashMap<usize, usize>,

    // interface-domain mapping
    pub node_itf_dom1_id: Vec<usize>,
    pub node_itf_dom2_id: Vec<usize>,
    pub node_dom1_itf_id: HashMap<usize, usize>,
    pub node_dom2_itf_id: HashMap<usize, usize>,
    pub elem_itf_dom1_id: Vec<usize>,
    pub elem_itf_dom2_id: Vec<usize>,
}

impl Interface {
    pub fn new(
        itf_id: usize,
        mesh: &Mesh,
        dom1: &Domain,
        dom2: &Domain,
        reg_id: usize,
    ) -> Result<Interface, FEChemError> {
        // initialize interface
        let mut itf = Interface::default();
        itf.itf_id = itf_id;
        itf.dom1_id = dom1.dom_id;
        itf.dom2_id = dom2.dom_id;

        // step 1: build interface-mesh mappings

        // create interface-mesh element mapping
        itf.elem_itf_mesh_id = mesh.reg1d_elem_id[reg_id].clone();
        for (elem_iid, &elem_mid) in itf.elem_itf_mesh_id.iter().enumerate() {
            itf.elem_mesh_itf_id.insert(elem_mid, elem_iid);
        }

        // build list of global node ids
        let mut node_set = HashSet::new();
        for &elem_mid in itf.elem_itf_mesh_id.iter() {
            for &node_mid in &mesh.elm1d_node_id[elem_mid] {
                node_set.insert(node_mid);
            }
        }

        // create global-local node mapping
        itf.node_itf_mesh_id = node_set.into_iter().collect::<Vec<usize>>();
        for (node_iid, &node_mid) in itf.node_itf_mesh_id.iter().enumerate() {
            itf.node_mesh_itf_id.insert(node_mid, node_iid);
        }

        // step 2: get data from mesh

        // get node coordinates
        for &node_mid in itf.node_itf_mesh_id.iter() {
            itf.node_x.push(mesh.node_x[node_mid]);
            itf.node_y.push(mesh.node_y[node_mid]);
        }
        itf.num_node = itf.node_itf_mesh_id.len();

        // get element data
        for &elem_mid in itf.elem_itf_mesh_id.iter() {
            // element type - copy from mesh
            itf.elem_node_num.push(mesh.elm1d_node_num[elem_mid]);

            // node ids - orient to match parent domain element edge traversal
            let mesh_nid = &mesh.elm1d_node_id[elem_mid];
            let d1_0 = dom1.node_mesh_dom_id[&mesh_nid[0]];
            let d1_1 = dom1.node_mesh_dom_id[&mesh_nid[1]];
            let (dom1_eid, dom1_n0, dom1_n1) = find_domain_edge(dom1, d1_0, d1_1)?;
            let m1_0 = dom1.node_dom_mesh_id[dom1_n0];
            let m1_1 = dom1.node_dom_mesh_id[dom1_n1];
            itf.elem_node1_id.push(vec![
                itf.node_mesh_itf_id[&m1_0],
                itf.node_mesh_itf_id[&m1_1],
            ]);

            let d2_0 = dom2.node_mesh_dom_id[&mesh_nid[0]];
            let d2_1 = dom2.node_mesh_dom_id[&mesh_nid[1]];
            let (dom2_eid, dom2_n0, dom2_n1) = find_domain_edge(dom2, d2_0, d2_1)?;
            let m2_0 = dom2.node_dom_mesh_id[dom2_n0];
            let m2_1 = dom2.node_dom_mesh_id[dom2_n1];
            itf.elem_node2_id.push(vec![
                itf.node_mesh_itf_id[&m2_0],
                itf.node_mesh_itf_id[&m2_1],
            ]);

            itf.elem_itf_dom1_id.push(dom1_eid);
            itf.elem_itf_dom2_id.push(dom2_eid);
        }
        itf.num_elem = itf.elem_itf_mesh_id.len();

        // step 3: build interface-domain mappings

        // iterate through interface nodes
        for (node_iid, &node_mid) in itf.node_itf_mesh_id.iter().enumerate() {
            let node_d1id = dom1.node_mesh_dom_id[&node_mid];
            itf.node_itf_dom1_id.push(node_d1id);
            itf.node_dom1_itf_id.insert(node_d1id, node_iid);

            let node_d2id = dom2.node_mesh_dom_id[&node_mid];
            itf.node_itf_dom2_id.push(node_d2id);
            itf.node_dom2_itf_id.insert(node_d2id, node_iid);
        }

        // result
        Ok(itf)
    }
}

/// Find a domain element edge matching `d0` and `d1`.
/// Returns the parent element id and node ids in CCW edge order.
fn find_domain_edge(
    dom: &Domain,
    d0: usize,
    d1: usize,
) -> Result<(usize, usize, usize), FEChemError> {
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
