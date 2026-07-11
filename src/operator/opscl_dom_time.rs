use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub wgt_id: usize, // mass scalar
    pub unk_id: usize, // unknown scalar
}

impl OpSclDomTime {
    pub fn new(dom_id: usize, wgt_id: usize, unk_id: usize) -> OpSclDomTime {
        // adds d(wgt * unk)/dt to LHS using backward Euler
        // in the apply function, factor is dt

        // create struct
        let mut oper_time = OpSclDomTime::default();
        oper_time.dom_id = dom_id;
        oper_time.wgt_id = wgt_id;
        oper_time.unk_id = unk_id;

        // result
        oper_time
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t_next: f64, dt: f64, factor: f64) {
        // backward Euler discretizes d(wgt * unk)/dt as (wgt * unk - wgt_prev * unk_prev) / dt
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the local matrix (wgt * unk)/dt and vector (wgt_prev * unk_prev)/dt are negative

        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let wgt = &vars.scl_dom[self.wgt_id];
        let unk = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local mass matrix

            // initialize local matrix
            let num_node = dom.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];
            let mut b_loc = vec![vec![0.0; num_node]; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];

            // assemble local mass matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        let n = tri3_eval(A_TRI3[qid], B_TRI3[qid]);
                        let t_prev = t_next - dt;
                        let wgt_val = wgt.compute_quad(vars, eid, qid, t_next);
                        let wgt_val_prev = wgt.compute_quad_prev(vars, eid, qid, t_prev);
                        let coeff = -factor * W_TRI3[qid] * wgt_val * jac_det[qid] / dt;
                        let coeff_prev = -factor * W_TRI3[qid] * wgt_val_prev * jac_det[qid] / dt;
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff * n[v] * n[j];
                                b_loc[v][j] += coeff_prev * n[v] * n[j];
                            }
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        let n = quad4_eval(A_QUAD4[qid], B_QUAD4[qid]);
                        let t_prev = t_next - dt;
                        let wgt_val = wgt.compute_quad(vars, eid, qid, t_next);
                        let wgt_val_prev = wgt.compute_quad_prev(vars, eid, qid, t_prev);
                        let coeff = -factor * W_QUAD4[qid] * wgt_val * jac_det[qid] / dt;
                        let coeff_prev = -factor * W_QUAD4[qid] * wgt_val_prev * jac_det[qid] / dt;
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff * n[v] * n[j];
                                b_loc[v][j] += coeff_prev * n[v] * n[j];
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
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk.node_dir[nid_v] {
                    continue;
                }

                // add current time step
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, self.unk_id, nid_j, a_loc[v][j]);
                }

                // add previous time step values
                let mut b_contrib = 0.0;
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    b_contrib += b_loc[v][j] * unk.node_prev[nid_j];
                }
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, b_contrib);

            }
        }
    }
}

impl OperatorBase for OpSclDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, _factor: f64) {
        // do not use this. use apply_time instead
        panic!("Used apply for OperatorTime. Must use apply_time instead.");
    }
}
