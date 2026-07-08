use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDivergence {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub vel_id: usize, // velocity
    pub unk_id: usize, // unknown scalar (pressure)
}

impl OpSclDomDivergence {
    pub fn new(dom_id: usize, den_id: usize, vel_id: usize, unk_id: usize) -> OpSclDomDivergence {
        // adds -div(den * vel) to RHS
        // LHS is 0 for steady state or d(den)/dt for transient

        // create struct
        let mut oper_diff = OpSclDomDivergence::default();
        oper_diff.dom_id = dom_id;
        oper_diff.den_id = den_id;
        oper_diff.vel_id = vel_id;
        oper_diff.unk_id = unk_id;

        // result
        oper_diff
    }
}

impl OperatorBase for OpSclDomDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // -(div(den * vel), test)_dom = -(den * vel * norm, test)_bnd + (den * vel, grad(test))_dom
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the sign of the local matrix entries is positive
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let den = &vars.scl_dom[self.den_id];
        let unk = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node_num[eid];
            let mut ax_loc = vec![vec![0.0; num_node]; num_node];  // x-component of velocity
            let mut ay_loc = vec![vec![0.0; num_node]; num_node];  // y-component of velocity

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];
            let gradn_x = &itg.gradn_x[eid];
            let gradn_y = &itg.gradn_y[eid];

            // assemble local matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        // get values
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let coeff = factor * W_TRI3[qid] * den_val * jac_det[qid];

                        // get test function values
                        let a = A_TRI3[qid];
                        let b = B_TRI3[qid];
                        let n = tri3_eval(a, b);

                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                ax_loc[v][j] += coeff * gradn_x[qid][v] * n[j];
                                ay_loc[v][j] += coeff * gradn_y[qid][v] * n[j];
                            }
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        // get values
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let coeff = factor * W_QUAD4[qid] * den_val * jac_det[qid];

                        // get test function values
                        let a = A_QUAD4[qid];
                        let b = B_QUAD4[qid];
                        let n = quad4_eval(a, b);
                        
                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                ax_loc[v][j] += coeff * gradn_x[qid][v] * n[j];
                                ay_loc[v][j] += coeff * gradn_y[qid][v] * n[j];
                            }
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
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 0, nid_j, ax_loc[v][j]);
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 1, nid_j, ay_loc[v][j]);
                }
            }
        }
    }
}
