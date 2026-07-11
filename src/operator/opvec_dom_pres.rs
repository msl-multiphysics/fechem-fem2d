use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomPressure {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize, // unknown vector
    pub pres_id: usize, // pressure
}

impl OpVecDomPressure {
    pub fn new(dom_id: usize, unk_id: usize, pres_id: usize) -> OpVecDomPressure {
        // adds -grad(pres) to RHS
        // LHS is 0 for steady state or d(unk)/dt for transient

        // create struct
        let mut oper_adv = OpVecDomPressure::default();
        oper_adv.dom_id = dom_id;
        oper_adv.pres_id = pres_id;
        oper_adv.unk_id = unk_id;

        // result
        oper_adv
    }
}

impl OperatorBase for OpVecDomPressure {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the sign of the local matrix entries is negative
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let unk = &vars.vec_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut ax_loc = vec![vec![0.0; num_node]; num_node];  // x momentum
            let mut ay_loc = vec![vec![0.0; num_node]; num_node];  // y momentum

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];
            let gradn_x = &itg.quad_gnx[eid];
            let gradn_y = &itg.quad_gny[eid];

            // assemble local matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        // get values
                        let coeff = -factor * W_TRI3[qid] *jac_det[qid];

                        // get test function values
                        let a = A_TRI3[qid];
                        let b = B_TRI3[qid];
                        let n = tri3_eval(a, b);

                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                ax_loc[v][j] += coeff * gradn_x[qid][j] * n[v];
                                ay_loc[v][j] += coeff * gradn_y[qid][j] * n[v];
                            }
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        // get values
                        let coeff = -factor * W_QUAD4[qid] * jac_det[qid];

                        // get test function values
                        let a = A_QUAD4[qid];
                        let b = B_QUAD4[qid];
                        let n = quad4_eval(a, b);
                        
                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                ax_loc[v][j] += coeff * gradn_x[qid][j] * n[v];
                                ay_loc[v][j] += coeff * gradn_y[qid][j] * n[v];
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
                    self.add_a_vecdom_scldom(vars, a_triplet, self.unk_id, 0, nid_v, self.pres_id, nid_j, ax_loc[v][j]);
                    self.add_a_vecdom_scldom(vars, a_triplet, self.unk_id, 1, nid_v, self.pres_id, nid_j, ay_loc[v][j]);
                }
            }
        }
    }
}
