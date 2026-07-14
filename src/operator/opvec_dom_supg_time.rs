use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomSupgTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub visc_id: usize, // viscosity

    // vectors
    pub unk_id: usize, // unknown vector
    pub drv_id: usize, // driving vector
}

impl OpVecDomSupgTime {
    pub fn new(dom_id: usize, den_id: usize, visc_id: usize, unk_id: usize, drv_id: usize) -> OpVecDomSupgTime {
        // adds the time-derivative SUPG residual contribution to momentum
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        //
        // den - density (den_i)
        // visc - viscosity (mu)
        // unk - unknown vector (v_i)
        // drv - driving vector (v_j)

        // create struct
        let mut oper_supg = OpVecDomSupgTime::default();
        oper_supg.dom_id = dom_id;
        oper_supg.den_id = den_id;
        oper_supg.visc_id = visc_id;
        oper_supg.unk_id = unk_id;
        oper_supg.drv_id = drv_id;

        // result
        oper_supg
    }

    fn compute_tau(&self, den_val: f64, visc_val: f64, drv_x: f64, drv_y: f64, jac_met: &[[f64; 2]; 2]) -> f64 {
        // metric-based steady SUPG time scale (same tau as OpVecDomSupgSteady)
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

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t_next: f64, dt: f64, factor: f64) {
        // time derivative is discretized using backward Euler
        // d(den * v)/dt = (den_next * v_next - den_curr * v_curr) / dt
        //
        // apply the weak form of the time-derivative SUPG term
        // tau * (v_j . grad(w), d(den * v)/dt)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add tau * (v_j . grad(w), (den_next * v_next)/dt)_dom -> add with negative sign to A
        // add -tau * (v_j . grad(w), (den_curr * v_curr)/dt)_dom to b

        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let visc_scl = &vars.scl_dom[self.visc_id];
        let unk_vec = &vars.vec_dom[self.unk_id];
        let drv_vec = &vars.vec_dom[self.drv_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix and vector

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

            // assemble local matrix
            for qid in 0..num_quad {
                // next time step (tau and streamline weight from current iterate)
                let den = den_scl.compute_quad(vars, eid, qid, t_next);
                let visc = visc_scl.compute_quad(vars, eid, qid, t_next);
                let (drv_x, drv_y) = drv_vec.compute_quad(vars, eid, qid, t_next);  // lag the driving vector by 1 iteration
                let tau = self.compute_tau(den, visc, drv_x, drv_y, &jac_met[qid]);
                let coeff = -factor * quad_w[qid] * tau * den * jac_det[qid] / dt;

                // current time step
                let t_curr = t_next - dt;
                let den_curr = den_scl.compute_quad_prev(vars, eid, qid, t_curr);
                let coeff_curr = -factor * quad_w[qid] * tau * den_curr * jac_det[qid] / dt;

                // load entries
                for v in 0..num_node {
                    let drv_grad_v = drv_x * quad_gnx[qid][v] + drv_y * quad_gny[qid][v];
                    for j in 0..num_node {
                        an_loc[v][j] += coeff * drv_grad_v * quad_n[qid][j];
                        ac_loc[v][j] += coeff_curr * drv_grad_v * quad_n[qid][j];
                    }
                }
            }

            // step 2: add to global matrix and vector

            // iterate over local matrix entries
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk_vec.node_dir[nid_v] {
                    continue;
                }

                // add next time step
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.unk_id, 0, nid_j, an_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.unk_id, 1, nid_j, an_loc[v][j]);
                }

                // add current time step
                let mut ac_sum_x = 0.0;
                let mut ac_sum_y = 0.0;
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    ac_sum_x += ac_loc[v][j] * unk_vec.node_value_prev_x[nid_j];
                    ac_sum_y += ac_loc[v][j] * unk_vec.node_value_prev_y[nid_j];
                }
                self.add_b_vecdom(vars, b_vec, self.unk_id, 0, nid_v, ac_sum_x);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 1, nid_v, ac_sum_y);
            }
        }
    }
}

impl OperatorBase for OpVecDomSupgTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, _factor: f64) {
        // do not use this. use apply_time instead
        panic!("Used apply for OpVecDomSupgTime. Must use apply_time instead.");
    }
}
