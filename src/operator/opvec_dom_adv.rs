use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomAdvection {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub unk_id: usize, // unknown vector
    pub drv_id: usize, // driving vector
}

impl OpVecDomAdvection {
    pub fn new(dom_id: usize, den_id: usize, unk_id: usize, drv_id: usize) -> OpVecDomAdvection {
        // adds -den * drv * grad(unk) to RHS
        // LHS is 0 for steady state or d(unk)/dt for transient

        // create struct
        let mut oper_adv = OpVecDomAdvection::default();
        oper_adv.dom_id = dom_id;
        oper_adv.den_id = den_id;
        oper_adv.drv_id = drv_id;
        oper_adv.unk_id = unk_id;

        // result
        oper_adv
    }
}

impl OperatorBase for OpVecDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the sign of the local matrix entries is negative
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let den = &vars.scl_dom[self.den_id];
        let unk = &vars.vec_dom[self.unk_id];
        let drv = &vars.vec_dom[self.drv_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];  // both x and y momentum have the same local matrix

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
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let (drv_x, drv_y) = drv.compute_quad(vars, eid, qid, t);  // lag the driving vector by 1 iteration
                        let coeff = -factor * W_TRI3[qid] * den_val * jac_det[qid];

                        // get test function values
                        let a = A_TRI3[qid];
                        let b = B_TRI3[qid];
                        let n = tri3_eval(a, b);

                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff * (drv_x * gradn_x[qid][j] + drv_y * gradn_y[qid][j]) * n[v];
                            }
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        // get values
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let (drv_x, drv_y) = drv.compute_quad(vars, eid, qid, t);  // lag the driving vector by 1 iteration
                        let coeff = -factor * W_QUAD4[qid] * den_val * jac_det[qid];

                        // get test function values
                        let a = A_QUAD4[qid];
                        let b = B_QUAD4[qid];
                        let n = quad4_eval(a, b);
                        
                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff * (drv_x * gradn_x[qid][j] + drv_y * gradn_y[qid][j]) * n[v];
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
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.unk_id, 0, nid_j, a_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.unk_id, 1, nid_j, a_loc[v][j]);
                }
            }
        }
    }
}
