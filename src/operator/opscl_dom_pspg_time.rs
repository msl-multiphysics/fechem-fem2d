use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomPspgTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub visc_id: usize, // viscosity
    pub unk_id: usize, // unknown scalar (pressure)

    // vectors
    pub vel_id: usize, // velocity
}

impl OpSclDomPspgTime {
    pub fn new(dom_id: usize, den_id: usize, visc_id: usize, vel_id: usize, unk_id: usize) -> OpSclDomPspgTime {
        // adds the time-derivative PSPG residual contribution to continuity
        // d(rho)/dt = -div(rho * v)
        //
        // den - density (rho)
        // visc - viscosity (mu)
        // vel - velocity (v)
        // unk - unknown scalar (equation added to rows of this scalar; e.g., pressure)

        // create struct
        let mut oper_pspg = OpSclDomPspgTime::default();
        oper_pspg.dom_id = dom_id;
        oper_pspg.den_id = den_id;
        oper_pspg.visc_id = visc_id;
        oper_pspg.vel_id = vel_id;
        oper_pspg.unk_id = unk_id;

        // result
        oper_pspg
    }

    fn compute_tau(&self, den_val: f64, visc_val: f64, vel_x: f64, vel_y: f64, jac_met: &[[f64; 2]; 2]) -> f64 {
        // metric-based steady stabilization time scale (same tau as OpSclDomPspgSteady)
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

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t_next: f64, dt: f64, factor: f64) {
        // time derivative is discretized using backward Euler
        // d(den * v)/dt = (den_next * v_next - den_curr * v_curr) / dt
        //
        // apply the weak form of the time-derivative PSPG term
        // tau * (grad(w), d(den * v)/dt)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add tau * (grad(w), (den_next * v_next)/dt)_dom -> add with negative sign to A
        // add -tau * (grad(w), (den_curr * v_curr)/dt)_dom to b

        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let visc_scl = &vars.scl_dom[self.visc_id];
        let vel_vec = &vars.vec_dom[self.vel_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix and vector

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut anx_loc = vec![vec![0.0; num_node]; num_node];  // next time step, x-velocity
            let mut any_loc = vec![vec![0.0; num_node]; num_node];  // next time step, y-velocity
            let mut acx_loc = vec![vec![0.0; num_node]; num_node];  // current time step, x-velocity
            let mut acy_loc = vec![vec![0.0; num_node]; num_node];  // current time step, y-velocity

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];
            let jac_met = &itgdom.jac_met[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                // next time step (tau from current iterate)
                let den = den_scl.compute_quad(vars, eid, qid, t_next);
                let visc = visc_scl.compute_quad(vars, eid, qid, t_next);
                let (vel_x, vel_y) = vel_vec.compute_quad(vars, eid, qid, t_next);  // lag the velocity by 1 iteration
                let tau = self.compute_tau(den, visc, vel_x, vel_y, &jac_met[qid]);
                let coeff = -factor * quad_w[qid] * tau * den * jac_det[qid] / dt;

                // current time step
                let t_curr = t_next - dt;
                let den_curr = den_scl.compute_quad_prev(vars, eid, qid, t_curr);
                let coeff_curr = -factor * quad_w[qid] * tau * den_curr * jac_det[qid] / dt;

                // load entries
                for v in 0..num_node {
                    for j in 0..num_node {
                        anx_loc[v][j] += coeff * quad_gnx[qid][v] * quad_n[qid][j];
                        any_loc[v][j] += coeff * quad_gny[qid][v] * quad_n[qid][j];
                        acx_loc[v][j] += coeff_curr * quad_gnx[qid][v] * quad_n[qid][j];
                        acy_loc[v][j] += coeff_curr * quad_gny[qid][v] * quad_n[qid][j];
                    }
                }
            }

            // step 2: add to global matrix and vector

            // iterate over local matrix entries
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
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 0, nid_j, anx_loc[v][j]);
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 1, nid_j, any_loc[v][j]);
                }

                // add current time step
                let mut ac_sum = 0.0;
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    ac_sum += acx_loc[v][j] * vel_vec.node_value_prev_x[nid_j];
                    ac_sum += acy_loc[v][j] * vel_vec.node_value_prev_y[nid_j];
                }
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, ac_sum);
            }
        }
    }
}

impl OperatorBase for OpSclDomPspgTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, _factor: f64) {
        // do not use this. use apply_time instead
        panic!("Used apply for OpSclDomPspgTime. Must use apply_time instead.");
    }
}
