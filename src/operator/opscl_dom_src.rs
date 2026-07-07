use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OperatorSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub src_id: usize, // source
    pub unk_id: usize, // unknown scalar
}

impl OperatorSource {
    pub fn new(dom_id: usize, src_id: usize, unk_id: usize) -> OperatorSource {
        // adds +src to RHS
        // LHS is 0 for steady state or d(unk)/dt for transient

        // create struct
        let mut oper_src = OperatorSource::default();
        oper_src.dom_id = dom_id;
        oper_src.src_id = src_id;
        oper_src.unk_id = unk_id;

        // result
        oper_src
    }
}

impl OperatorBase for OperatorSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the sign of the local matrix entries is negative
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let src = &vars.scl_dom[self.src_id];
        let unk = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local vector

            // initialize local vectors
            let num_node = dom.elem_node_num[eid];
            let mut b_loc = vec![0.0; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];

            // assemble local vector
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        let n = tri3_eval(A_TRI3[qid], B_TRI3[qid]);
                        let src_val = src.compute_quad(vars, eid, qid, t);
                        let coeff = -factor * W_TRI3[qid] * src_val * jac_det[qid];
                        for v in 0..num_node {
                            b_loc[v] += coeff * n[v];
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        let n = quad4_eval(A_QUAD4[qid], B_QUAD4[qid]);
                        let src_val = src.compute_quad(vars, eid, qid, t);
                        let coeff = -factor * W_QUAD4[qid] * src_val * jac_det[qid];
                        for v in 0..num_node {
                            b_loc[v] += coeff * n[v];
                        }
                    }
                }
                _ => {
                    panic!("Invalid element type");
                }
            }

            // step 2: assemble global vector

            // iterate over local vector entries
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk.node_dir[nid_v] {
                    continue;
                }

                // add to global vector
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, b_loc[v]);
            }
        }
    }
}
