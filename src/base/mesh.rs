use crate::base::error::FEChemError;
use crate::base::read_gmsh::read_gmsh_mesh;

#[derive(Default)]
pub struct Mesh {
    // node data
    pub num_node: usize,
    pub node_x: Vec<f64>,
    pub node_y: Vec<f64>,

    // 2d element data
    pub num_elm2d: usize,               // number of 2d elements
    pub num_reg2d: usize,               // number of 2d regions
    pub elm2d_node_num: Vec<usize>,     // number of nodes (3 - tri; 4 - quad)
    pub elm2d_node_id: Vec<Vec<usize>>, // node ids of elements
    pub reg2d_elem_id: Vec<Vec<usize>>, // reg2d_elem_id[reg_id] -> list of element ids belonging to region

    // 1d element data
    pub num_elm1d: usize,
    pub num_reg1d: usize,           // number of 1d regions
    pub elm1d_node_num: Vec<usize>, // number of nodes (2 - lin)
    pub elm1d_node_id: Vec<Vec<usize>>,
    pub reg1d_elem_id: Vec<Vec<usize>>, // reg1d_elem_id[reg_id] -> list of element ids belonging to region
}

impl Mesh {
    pub fn new(file_path: String) -> Result<Mesh, FEChemError> {
        read_gmsh_mesh(&file_path)
    }

    pub fn new_from_bounds(
        x_min: f64,
        y_min: f64,
        x_max: f64,
        y_max: f64,
        num_elm_x: usize,
        num_elm_y: usize,
    ) -> Result<Mesh, FEChemError> {
        // initialize Mesh
        let mut dom = Mesh::default();

        // compute discretization
        let dx = (x_max - x_min) / (num_elm_x as f64);
        let dy = (y_max - y_min) / (num_elm_y as f64);
        let stride = num_elm_x + 1; // nodes per row in the structured grid

        // create nodes
        for j in 0..=num_elm_y {
            for i in 0..=num_elm_x {
                let x = x_min + i as f64 * dx;
                let y = y_min + j as f64 * dy;
                dom.node_x.push(x);
                dom.node_y.push(y);
            }
        }
        dom.num_node = stride * (num_elm_y + 1);

        // create 2d elements
        for j in 0..num_elm_y {
            for i in 0..num_elm_x {
                let nid0 = i + j * stride; // node id 0 - bottom left
                let nid1 = nid0 + 1; // node id 1 - bottom right
                let nid2 = nid0 + stride + 1; // node id 2 - top right
                let nid3 = nid0 + stride; // node id 3 - top left
                dom.elm2d_node_num.push(4); // quad4
                dom.elm2d_node_id.push(vec![nid0, nid1, nid2, nid3]);
            }
        }
        dom.num_elm2d = num_elm_x * num_elm_y;
        dom.num_reg2d = 1;
        dom.reg2d_elem_id.push((0..num_elm_x * num_elm_y).collect());

        // create 1d elements
        // same direction as 2d elements
        for j in (0..num_elm_y).rev() {
            // left
            let nid0 = (j + 1) * stride; // node id 0 - top
            let nid1 = j * stride; // node id 1 - bottom
            dom.elm1d_node_num.push(2); // lin2
            dom.elm1d_node_id.push(vec![nid0, nid1]);
        }
        for j in 0..num_elm_y {
            // right
            let nid0 = num_elm_x + j * stride; // node id 0 - bottom
            let nid1 = num_elm_x + (j + 1) * stride; // node id 1 - top
            dom.elm1d_node_num.push(2); // lin2
            dom.elm1d_node_id.push(vec![nid0, nid1]);
        }
        for i in 0..num_elm_x {
            // bottom
            let nid0 = i; // node id 0 - left
            let nid1 = i + 1; // node id 1 - right
            dom.elm1d_node_num.push(2); // lin2
            dom.elm1d_node_id.push(vec![nid0, nid1]);
        }
        for i in (0..num_elm_x).rev() {
            // top
            let nid0 = num_elm_y * stride + i + 1; // node id 0 - right
            let nid1 = num_elm_y * stride + i; // node id 1 - left
            dom.elm1d_node_num.push(2); // lin2
            dom.elm1d_node_id.push(vec![nid0, nid1]);
        }
        dom.num_elm1d = 2 * (num_elm_x + num_elm_y);
        dom.num_reg1d = 4;
        dom.reg1d_elem_id.push((0..num_elm_y).collect());
        dom.reg1d_elem_id.push((num_elm_y..2 * num_elm_y).collect());
        dom.reg1d_elem_id
            .push((2 * num_elm_y..2 * num_elm_y + num_elm_x).collect());
        dom.reg1d_elem_id
            .push((2 * num_elm_y + num_elm_x..2 * num_elm_y + 2 * num_elm_x).collect());

        // result
        Ok(dom)
    }
}
