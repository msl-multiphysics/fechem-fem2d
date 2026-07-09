use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomSupg {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub visc_id: usize, // viscosity
    pub pres_id: usize, // pressure
    pub fce_id: usize, // body force

    // vectors
    pub unk_id: usize, // unknown vector
    pub drv_id: usize, // driving vector
}

impl OpVecDomSupg {
    pub fn new(dom_id: usize, den_id: usize, visc_id: usize, unk_id: usize, pres_id: usize, fce_id: usize, drv_id: usize) -> OpVecDomSupg {
        // adds SUPG stabilization for the momentum equation to RHS
        // LHS is 0 for steady state or d(unk)/dt for transient

        // create struct
        let mut oper_supg = OpVecDomSupg::default();
        oper_supg.dom_id = dom_id;
        oper_supg.den_id = den_id;
        oper_supg.visc_id = visc_id;
        oper_supg.unk_id = unk_id;
        oper_supg.pres_id = pres_id;
        oper_supg.fce_id = fce_id;
        oper_supg.drv_id = drv_id;

        // result
        oper_supg
    }

    fn compute_tau(&self, den_val: f64, visc_val: f64, drv_x: f64, drv_y: f64, jac_met: &[[f64; 2]; 2]) -> f64 {
        // metric-based steady SUPG time scale
        let g00 = jac_met[0][0];
        let g01 = jac_met[0][1];
        let g10 = jac_met[1][0];
        let g11 = jac_met[1][1];
        let adv = (drv_x * (g00 * drv_x + g01 * drv_y) + drv_y * (g10 * drv_x + g11 * drv_y)).max(0.0).sqrt();
        let diff = if den_val.abs() > 1e-30 {
            (visc_val / den_val).abs() * (g00 * g00 + g01 * g01 + g10 * g10 + g11 * g11).sqrt()
        } else {
            0.0
        };

        1.0 / (2.0 * adv + 4.0 * diff + 1e-30)
    }
}

impl OperatorBase for OpVecDomSupg {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // assume that A (in Ax = b) is the RHS of the PDE; b is on the LHS of the PDE
        // therefore, the sign of the local matrix entries is negative
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let den = &vars.scl_dom[self.den_id];
        let visc = &vars.scl_dom[self.visc_id];
        let fce = &vars.vec_dom[self.fce_id];
        let unk = &vars.vec_dom[self.unk_id];
        let drv = &vars.vec_dom[self.drv_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix and vector

            // initialize local matrices
            let num_node = dom.elem_node_num[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];  // both x and y momentum have the same local matrix
            let mut axp_loc = vec![vec![0.0; num_node]; num_node];  // pressure coupling in x momentum
            let mut ayp_loc = vec![vec![0.0; num_node]; num_node];  // pressure coupling in y momentum
            let mut bx_loc = vec![0.0; num_node];
            let mut by_loc = vec![0.0; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];
            let jac_met = &itg.jac_met[eid];
            let gradn_x = &itg.gradn_x[eid];
            let gradn_y = &itg.gradn_y[eid];

            // assemble local matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        // get values
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let visc_val = visc.compute_quad(vars, eid, qid, t);
                        let (drv_x, drv_y) = drv.compute_quad(vars, eid, qid, t);  // lag the driving vector by 1 iteration
                        let (fce_x, fce_y) = fce.compute_quad(vars, eid, qid, t);
                        let tau = self.compute_tau(den_val, visc_val, drv_x, drv_y, &jac_met[qid]);
                        let coeff = -factor * W_TRI3[qid] * tau * jac_det[qid];

                        // add to local matrix
                        for v in 0..num_node {
                            let drv_grad_v = drv_x * gradn_x[qid][v] + drv_y * gradn_y[qid][v];
                            for j in 0..num_node {
                                let drv_grad_j = drv_x * gradn_x[qid][j] + drv_y * gradn_y[qid][j];
                                a_loc[v][j] += coeff * drv_grad_v * den_val * drv_grad_j;
                                axp_loc[v][j] += coeff * drv_grad_v * gradn_x[qid][j];
                                ayp_loc[v][j] += coeff * drv_grad_v * gradn_y[qid][j];
                            }
                            bx_loc[v] += coeff * drv_grad_v * fce_x;
                            by_loc[v] += coeff * drv_grad_v * fce_y;
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        // get values
                        let den_val = den.compute_quad(vars, eid, qid, t);
                        let visc_val = visc.compute_quad(vars, eid, qid, t);
                        let (drv_x, drv_y) = drv.compute_quad(vars, eid, qid, t);  // lag the driving vector by 1 iteration
                        let (fce_x, fce_y) = fce.compute_quad(vars, eid, qid, t);
                        let tau = self.compute_tau(den_val, visc_val, drv_x, drv_y, &jac_met[qid]);
                        let coeff = -factor * W_QUAD4[qid] * tau * jac_det[qid];
                        
                        // add to local matrix
                        for v in 0..num_node {
                            let drv_grad_v = drv_x * gradn_x[qid][v] + drv_y * gradn_y[qid][v];
                            for j in 0..num_node {
                                let drv_grad_j = drv_x * gradn_x[qid][j] + drv_y * gradn_y[qid][j];
                                a_loc[v][j] += coeff * drv_grad_v * den_val * drv_grad_j;
                                axp_loc[v][j] += coeff * drv_grad_v * gradn_x[qid][j];
                                ayp_loc[v][j] += coeff * drv_grad_v * gradn_y[qid][j];
                            }
                            bx_loc[v] += coeff * drv_grad_v * fce_x;
                            by_loc[v] += coeff * drv_grad_v * fce_y;
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
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.unk_id, 0, nid_j, a_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.unk_id, 1, nid_j, a_loc[v][j]);
                    self.add_a_vecdom_scldom(vars, a_triplet, self.unk_id, 0, nid_v, self.pres_id, nid_j, axp_loc[v][j]);
                    self.add_a_vecdom_scldom(vars, a_triplet, self.unk_id, 1, nid_v, self.pres_id, nid_j, ayp_loc[v][j]);
                }
                self.add_b_vecdom(vars, b_vec, self.unk_id, 0, nid_v, bx_loc[v]);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 1, nid_v, by_loc[v]);
            }
        }
    }
}
