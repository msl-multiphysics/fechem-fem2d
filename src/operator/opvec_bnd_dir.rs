use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndDirichlet {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub val_id: usize, // prescribed value
    pub unk_id: usize, // unknown vector
}

impl OpVecBndDirichlet {
    pub fn new(bnd_id: usize, val_id: usize, unk_id: usize) -> OpVecBndDirichlet {
        // imposes unk = val on the unknown scalar

        // create struct
        let mut oper_dir = OpVecBndDirichlet::default();
        oper_dir.bnd_id = bnd_id;
        oper_dir.val_id = val_id;
        oper_dir.unk_id = unk_id;

        // result
        oper_dir
    }
}

impl OperatorBase for OpVecBndDirichlet {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let val = &vars.vec_bnd[self.val_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // step 1: get local nodes
            let num_node = bnd.elem_node[eid];
            let node_id = &bnd.elem_node_id[eid];

            // step 2: assemble global matrix and vector

            // iterate over element nodes
            for v in 0..num_node {
                // get node ids
                let nid_bnd = node_id[v];
                let nid_dom = bnd.node_bnd_dom_id[nid_bnd];

                // get properties
                let val_x = val.node_value_x[nid_bnd];
                let val_y = val.node_value_y[nid_bnd];

                // impose dirichlet BC
                self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_dom, self.unk_id, 0, nid_dom, 1.0);
                self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_dom, self.unk_id, 1, nid_dom, 1.0);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 0, nid_dom, factor * val_x);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 1, nid_dom, factor * val_y);

            }
        }
    }
}
