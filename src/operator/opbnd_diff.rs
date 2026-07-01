use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OperatorNeumannDiffusion {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub flx_id: usize,  // flux
    pub unk_id: usize,  // unknown scalar
}

impl OperatorNeumannDiffusion {
    pub fn new(bnd_id: usize, flx_id: usize, unk_id: usize) -> OperatorNeumannDiffusion {
        // applies flux to the unknown scalar

        // create struct
        let mut oper_diff = OperatorNeumannDiffusion::default();
        oper_diff.bnd_id = bnd_id;
        oper_diff.flx_id = flx_id;
        oper_diff.unk_id = unk_id;

        // result
        oper_diff
    }
}

impl OperatorBase for OperatorNeumannDiffusion {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, factor: f64) {    
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let itg = &vars.itg_bnd[self.bnd_id];
        let flx = &vars.scl_bnd[self.flx_id];
        let unk = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // step 1: assemble local vector

            // initialize local vectors
            let num_node = bnd.elem_node[eid];
            let mut b_loc = vec![0.0; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];

            // get properties
            let flx = &flx.quad_value[eid];

            // assemble local vector
            match num_node {
                2 => {
                    for qid in 0..num_quad {
                        let n = lin2_eval(A_LIN2[qid]);
                        let coeff = -factor * W_LIN2[qid] * jac_det[qid].sqrt();
                        for v in 0..num_node {
                            b_loc[v] += coeff * flx[qid] * n[v];
                        }
                    }
                }
                _ => {panic!("Invalid element type");}
            }

            // step 2: assemble global vector

            // iterate over local vector entries
            let node_id = &bnd.elem_node_id[eid];
            for v in 0..num_node {
                // get node ids
                let nid_bnd = node_id[v];
                let nid_dom = bnd.node_bnd_dom_id[nid_bnd];

                // skip if dirichlet BC
                if unk.node_dir[nid_dom] {
                    continue;
                }

                // add to global vector
                self.add_b(vars, b_vec, self.unk_id, nid_dom, b_loc[v]);
            }

        }
    
    }
}
