use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub src_id: usize, // source
    pub unk_id: usize, // unknown vector

}

impl OpVecDomSource {
    pub fn new(dom_id: usize, src_id: usize, unk_id: usize) -> OpVecDomSource {
        // adds +src to RHS
        // LHS is 0 for steady state or d(unk)/dt for transient

        // create struct
        let mut oper_adv = OpVecDomSource::default();
        oper_adv.dom_id = dom_id;
        oper_adv.src_id = src_id;
        oper_adv.unk_id = unk_id;

        // result
        oper_adv
    }
}

impl OperatorBase for OpVecDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the sign of the local matrix entries is negative
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let src = &vars.vec_dom[self.src_id];
        let unk = &vars.vec_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut bx_loc = vec![0.0; num_node];
            let mut by_loc = vec![0.0; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];

            // assemble local matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        // get values
                        let (src_x, src_y) = src.compute_quad(vars, eid, qid, t);
                        let coeff = -factor * W_TRI3[qid] * jac_det[qid];
                        
                        // get test function values
                        let a = A_TRI3[qid];
                        let b = B_TRI3[qid];
                        let n = tri3_eval(a, b);

                        // add to local matrix
                        for v in 0..num_node {
                            bx_loc[v] += coeff * src_x * n[v];
                            by_loc[v] += coeff * src_y * n[v];
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        // get values
                        let (src_x, src_y) = src.compute_quad(vars, eid, qid, t);
                        let coeff = -factor * W_QUAD4[qid] * jac_det[qid];

                        // get test function values
                        let a = A_QUAD4[qid];
                        let b = B_QUAD4[qid];
                        let n = quad4_eval(a, b);

                        // add to local matrix
                        for v in 0..num_node {
                            bx_loc[v] += coeff * src_x * n[v];
                            by_loc[v] += coeff * src_y * n[v];
                        }
                    }
                }
                _ => {
                    panic!("Invalid element type");
                }
            }

            // step 2: assemble global matrix

            // iterate over local matrix entries
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk.node_dir[nid_v] {
                    continue;
                }
                
                // add to global matrix
                self.add_b_vecdom(vars, b_vec, self.unk_id, 0, nid_v, bx_loc[v]);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 1, nid_v, by_loc[v]);
            }
        }
    }
}
