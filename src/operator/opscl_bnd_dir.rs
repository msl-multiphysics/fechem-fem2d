use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct OpSclBndDirichlet {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub val_id: usize, // prescribed value
    pub unk_id: usize, // unknown scalar
}

impl OpSclBndDirichlet {
    pub fn new(bnd_id: usize, val_id: usize, unk_id: usize) -> OpSclBndDirichlet {
        // applies a Dirichlet BC
        // c_i = u_i
        // 
        // val - prescribed value (u_i)
        // unk - unknown scalar (c_i)

        // create struct
        let mut oper_dir = OpSclBndDirichlet::default();
        oper_dir.bnd_id = bnd_id;
        oper_dir.val_id = val_id;
        oper_dir.unk_id = unk_id;

        // result
        oper_dir
    }

    pub fn apply_initial(&self, vars: &Variables, sum: &mut HashMap<(usize, usize), f64>, count: &mut HashMap<(usize, usize), usize>) {
        // accumulate prescribed values for initial-condition projection
        // visits the same (element, node) pairs as apply so that corner
        // conflicts and shared edge nodes average the same way as Ax = b
        // key: (unk_id, nid_dom)

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let val = &vars.scl_bnd[self.val_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // iterate over nodes in element
            let num_node = bnd.elem_node[eid];
            let node_id = &bnd.elem_node_id[eid];
            for v in 0..num_node {
                // get node ids
                let nid_bnd = node_id[v];
                let nid_dom = bnd.node_bnd_dom_id[nid_bnd];

                // accumulate prescribed value
                let key = (self.unk_id, nid_dom);
                *sum.entry(key).or_insert(0.0) += val.node_value[nid_bnd];
                *count.entry(key).or_insert(0) += 1;
            }
        }
    }
}

impl OperatorBase for OpSclBndDirichlet {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let val = &vars.scl_bnd[self.val_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // iterate over nodes in element
            let num_node = bnd.elem_node[eid];
            let node_id = &bnd.elem_node_id[eid];
            for v in 0..num_node {
                // get node ids
                let nid_bnd = node_id[v];
                let nid_dom = bnd.node_bnd_dom_id[nid_bnd];

                // get properties
                let val = val.node_value[nid_bnd];

                // impose dirichlet BC
                self.add_a_scldom(vars, a_triplet, self.unk_id, nid_dom, self.unk_id, nid_dom, 1.0);
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_dom, factor * val);
            }
        }

    }
}
