use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomPspg {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub visc_id: usize, // viscosity
    pub pres_id: usize, // pressure
    pub fce_id: usize, // body force
    pub unk_id: usize, // unknown scalar (pressure)

    // vectors
    pub vel_id: usize, // velocity
}

impl OpSclDomPspg {
    pub fn new(dom_id: usize, den_id: usize, visc_id: usize, vel_id: usize, pres_id: usize, fce_id: usize, unk_id: usize) -> OpSclDomPspg {
        // adds PSPG stabilization for the continuity equation to RHS
        // LHS is 0 for steady state or d(den)/dt for transient

        // create struct
        let mut oper_pspg = OpSclDomPspg::default();
        oper_pspg.dom_id = dom_id;
        oper_pspg.den_id = den_id;
        oper_pspg.visc_id = visc_id;
        oper_pspg.vel_id = vel_id;
        oper_pspg.pres_id = pres_id;
        oper_pspg.fce_id = fce_id;
        oper_pspg.unk_id = unk_id;

        // result
        oper_pspg
    }

    fn compute_tau(&self, den_val: f64, visc_val: f64, vel_x: f64, vel_y: f64, jac_met: &[[f64; 2]; 2]) -> f64 {
        // metric-based steady stabilization time scale
        let g00 = jac_met[0][0];
        let g01 = jac_met[0][1];
        let g10 = jac_met[1][0];
        let g11 = jac_met[1][1];
        let adv = (vel_x * (g00 * vel_x + g01 * vel_y) + vel_y * (g10 * vel_x + g11 * vel_y)).max(0.0).sqrt();
        let diff = if den_val.abs() > 1e-30 {
            (visc_val / den_val).abs() * (g00 * g00 + g01 * g01 + g10 * g10 + g11 * g11).sqrt()
        } else {
            0.0
        };

        1.0 / (2.0 * adv + 4.0 * diff + 1e-30)
    }
}

impl OperatorBase for OpSclDomPspg {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the sign of the local matrix entries is negative
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let den = &vars.scl_dom[self.den_id];
        let visc = &vars.scl_dom[self.visc_id];
        let fce = &vars.vec_dom[self.fce_id];
        let vel = &vars.vec_dom[self.vel_id];
        let unk = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix and vector

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut ax_loc = vec![vec![0.0; num_node]; num_node];  // x-component of velocity
            let mut ay_loc = vec![vec![0.0; num_node]; num_node];  // y-component of velocity
            let mut ap_loc = vec![vec![0.0; num_node]; num_node];  // pressure
            let mut b_loc = vec![0.0; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];
            let jac_met = &itg.jac_met[eid];
            let gradn_x = &itg.quad_gnx[eid];
            let gradn_y = &itg.quad_gny[eid];

            // assemble local matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        // get values
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let visc_val = visc.compute_quad(vars, eid, qid, t);
                        let (vel_x, vel_y) = vel.compute_quad(vars, eid, qid, t);  // lag the velocity by 1 iteration
                        let (fce_x, fce_y) = fce.compute_quad(vars, eid, qid, t);
                        let tau = self.compute_tau(den_val, visc_val, vel_x, vel_y, &jac_met[qid]);
                        let coeff = -factor * W_TRI3[qid] * tau * jac_det[qid];

                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                let vel_grad_j = vel_x * gradn_x[qid][j] + vel_y * gradn_y[qid][j];
                                ax_loc[v][j] += coeff * gradn_x[qid][v] * den_val * vel_grad_j;
                                ay_loc[v][j] += coeff * gradn_y[qid][v] * den_val * vel_grad_j;
                                ap_loc[v][j] += coeff * (gradn_x[qid][v] * gradn_x[qid][j] + gradn_y[qid][v] * gradn_y[qid][j]);
                            }
                            b_loc[v] += coeff * (gradn_x[qid][v] * fce_x + gradn_y[qid][v] * fce_y);
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        // get values
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let visc_val = visc.compute_quad(vars, eid, qid, t);
                        let (vel_x, vel_y) = vel.compute_quad(vars, eid, qid, t);  // lag the velocity by 1 iteration
                        let (fce_x, fce_y) = fce.compute_quad(vars, eid, qid, t);
                        let tau = self.compute_tau(den_val, visc_val, vel_x, vel_y, &jac_met[qid]);
                        let coeff = -factor * W_QUAD4[qid] * tau * jac_det[qid];
                        
                        // add to local matrix
                        for v in 0..num_node {
                            for j in 0..num_node {
                                let vel_grad_j = vel_x * gradn_x[qid][j] + vel_y * gradn_y[qid][j];
                                ax_loc[v][j] += coeff * gradn_x[qid][v] * den_val * vel_grad_j;
                                ay_loc[v][j] += coeff * gradn_y[qid][v] * den_val * vel_grad_j;
                                ap_loc[v][j] += coeff * (gradn_x[qid][v] * gradn_x[qid][j] + gradn_y[qid][v] * gradn_y[qid][j]);
                            }
                            b_loc[v] += coeff * (gradn_x[qid][v] * fce_x + gradn_y[qid][v] * fce_y);
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
                
                // add to global matrix and vector
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 0, nid_j, ax_loc[v][j]);
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 1, nid_j, ay_loc[v][j]);
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, self.pres_id, nid_j, ap_loc[v][j]);
                }
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, b_loc[v]);
            }
        }
    }
}
