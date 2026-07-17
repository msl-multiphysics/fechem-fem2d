use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomSupgTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub wgt_id: usize, // weighting scalar
    pub tau_diff_id: usize, // effective diffusion coefficient used in tau
    pub unk_id: usize, // unknown scalar

    // vectors
    pub vel_id: usize, // velocity vector
}

impl OpSclDomSupgTime {
    pub fn new(dom_id: usize, wgt_id: usize, tau_diff_id: usize, vel_id: usize, unk_id: usize) -> OpSclDomSupgTime {
        // adds the time-derivative SUPG residual contribution to scalar transport
        // d(m_i * c_i)/dt = -div(m_i * c_i * v - D_ij * grad(c_j)) + R_i
        //
        // wgt - weighting scalar (m_i)
        // tau_diff - effective diffusion coefficient used in tau (normally D_ii)
        // vel - velocity vector (v)
        // unk - unknown scalar (c_i)

        // create struct
        let mut oper_supg = OpSclDomSupgTime::default();
        oper_supg.dom_id = dom_id;
        oper_supg.wgt_id = wgt_id;
        oper_supg.tau_diff_id = tau_diff_id;
        oper_supg.vel_id = vel_id;
        oper_supg.unk_id = unk_id;

        // result
        oper_supg
    }

    fn compute_tau(&self, wgt_val: f64, diff_val: f64, vel_x: f64, vel_y: f64, jac_met: &[[f64; 2]; 2]) -> f64 {
        // metric-based steady SUPG time scale (same tau as OpSclDomSupgSteady)
        let g00 = jac_met[0][0];
        let g01 = jac_met[0][1];
        let g10 = jac_met[1][0];
        let g11 = jac_met[1][1];
        let adv = (vel_x * (g00 * vel_x + g01 * vel_y) + vel_y * (g10 * vel_x + g11 * vel_y)).max(0.0).sqrt();
        let diff = if wgt_val.abs() > 1e-30 {
            (diff_val / wgt_val).abs() * (g00 * g00 + g01 * g01 + g10 * g10 + g11 * g11).sqrt()
        } else {
            0.0
        };

        1.0 / (2.0 * adv + 4.0 * diff + 1e-30)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t_next: f64, dt: f64, factor: f64) {
        // time derivative is discretized using backward Euler
        // d(m * c)/dt = (m_next * c_next - m_curr * c_curr) / dt
        //
        // applies tau * (v . grad(w), d(m * c)/dt)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add the next-time contribution to A and the previous-time contribution to b

        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let wgt_scl = &vars.scl_dom[self.wgt_id];
        let tau_diff_scl = &vars.scl_dom[self.tau_diff_id];
        let unk_scl = &vars.scl_dom[self.unk_id];
        let vel_vec = &vars.vec_dom[self.vel_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrices

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut an_loc = vec![vec![0.0; num_node]; num_node];  // next time step
            let mut ac_loc = vec![vec![0.0; num_node]; num_node];  // current time step

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];
            let jac_met = &itgdom.jac_met[eid];

            // assemble local matrices
            for qid in 0..num_quad {
                // next time step (tau and streamline weight from current iterate)
                let wgt = wgt_scl.compute_quad(vars, eid, qid, t_next);
                let tau_diff = tau_diff_scl.compute_quad(vars, eid, qid, t_next);
                let (vel_x, vel_y) = vel_vec.compute_quad(vars, eid, qid, t_next);  // lag the velocity by 1 iteration
                let tau = self.compute_tau(wgt, tau_diff, vel_x, vel_y, &jac_met[qid]);
                let coeff = -factor * quad_w[qid] * tau * wgt * jac_det[qid] / dt;

                // current time step
                let t_curr = t_next - dt;
                let wgt_curr = wgt_scl.compute_quad_prev(vars, eid, qid, t_curr);
                let coeff_curr = -factor * quad_w[qid] * tau * wgt_curr * jac_det[qid] / dt;

                // load entries
                for v in 0..num_node {
                    let vel_grad_v = vel_x * quad_gnx[qid][v] + vel_y * quad_gny[qid][v];
                    for j in 0..num_node {
                        an_loc[v][j] += coeff * vel_grad_v * quad_n[qid][j];
                        ac_loc[v][j] += coeff_curr * vel_grad_v * quad_n[qid][j];
                    }
                }
            }

            // step 2: add to global matrix and vector
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk_scl.node_dir[nid_v] {
                    continue;
                }

                // add next time step
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, self.unk_id, nid_j, an_loc[v][j]);
                }

                // add current time step
                let mut ac_sum = 0.0;
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    ac_sum += ac_loc[v][j] * unk_scl.node_value_prev[nid_j];
                }
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, ac_sum);
            }
        }
    }
}

impl OperatorBase for OpSclDomSupgTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, _factor: f64) {
        // do not use this. use apply_time instead
        panic!("Used apply for OpSclDomSupgTime. Must use apply_time instead.");
    }
}
