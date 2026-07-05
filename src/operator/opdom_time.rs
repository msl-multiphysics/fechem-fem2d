use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OperatorTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize, // unknown scalar
}

impl OperatorTime {
    pub fn new(dom_id: usize, unk_id: usize) -> OperatorTime {
        // adds d(unk)/dt to LHS using backward Euler
        // in the apply function, factor is dt

        // create struct
        let mut oper_time = OperatorTime::default();
        oper_time.dom_id = dom_id;
        oper_time.unk_id = unk_id;

        // result
        oper_time
    }
}

impl OperatorBase for OperatorTime {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // factor is dt; backward Euler discretizes d(unk)/dt as (unk - unk_prev) / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let unk = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local mass matrix

            // initialize local matrix
            let num_node = dom.elem_node_num[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];

            // assemble local mass matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        let n = tri3_eval(A_TRI3[qid], B_TRI3[qid]);
                        let coeff = W_TRI3[qid] * jac_det[qid];
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff * n[v] * n[j];
                            }
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        let n = quad4_eval(A_QUAD4[qid], B_QUAD4[qid]);
                        let coeff = W_QUAD4[qid] * jac_det[qid];
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff * n[v] * n[j];
                            }
                        }
                    }
                }
                _ => {
                    panic!("Invalid element type");
                }
            }

            // step 2: assemble global matrix and vector

            // iterate over local matrix entries
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                // get node ids
                let nid_v = node_id[v];

                // skip if dirichlet BC
                if unk.node_dir[nid_v] {
                    continue;
                }

                // add previous time step values to RHS: (M / dt) * unk_prev
                let mut b_contrib = 0.0;
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    b_contrib += a_loc[v][j] * unk.node_prev[nid_j];
                }
                self.add_b_scl(vars, b_vec, self.unk_id, nid_v, b_contrib / factor);

                // add current time step unknowns to LHS: (M / dt) * unk
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_sclscl(vars, a_triplet, self.unk_id, nid_v, self.unk_id, nid_j, a_loc[v][j] / factor);
                }
            }
        }
    }
}
