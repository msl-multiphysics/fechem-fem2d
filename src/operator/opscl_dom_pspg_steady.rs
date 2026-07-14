use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomPspgSteady {
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

impl OpSclDomPspgSteady {
    pub fn new(dom_id: usize, den_id: usize, visc_id: usize, vel_id: usize, pres_id: usize, fce_id: usize, unk_id: usize) -> OpSclDomPspgSteady {
        // adds steady PSPG stabilization to the continuity equation
        // d(rho)/dt = -div(rho * v)
        // 
        // den - density (rho)
        // visc - viscosity (mu)
        // vel - velocity (v)
        // pres - pressure (p)
        // fce - body force (f)
        // unk - unknown scalar (equation added to rows of this scalar; e.g., pressure)

        // create struct
        let mut oper_pspg = OpSclDomPspgSteady::default();
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

impl OperatorBase for OpSclDomPspgSteady {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the steady PSPG stabilization term
        // tau * (grad(w), rho * v . grad(v) + grad(p) - f)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add the PSPG stabilization contributions to A and b
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let visc_scl = &vars.scl_dom[self.visc_id];
        let fce_vec = &vars.vec_dom[self.fce_id];
        let vel_vec = &vars.vec_dom[self.vel_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix and vector

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut ax_loc = vec![vec![0.0; num_node]; num_node];  // x-component of velocity
            let mut ay_loc = vec![vec![0.0; num_node]; num_node];  // y-component of velocity
            let mut ap_loc = vec![vec![0.0; num_node]; num_node];  // pressure
            let mut b_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];
            let jac_met = &itgdom.jac_met[eid];

            // assemble local matrix and vector
            for qid in 0..num_quad {
                let den = den_scl.compute_quad(vars, eid, qid, t);
                let visc = visc_scl.compute_quad(vars, eid, qid, t);
                let (vel_x, vel_y) = vel_vec.compute_quad(vars, eid, qid, t);  // lag the velocity by 1 iteration
                let (fce_x, fce_y) = fce_vec.compute_quad(vars, eid, qid, t);
                let tau = self.compute_tau(den, visc, vel_x, vel_y, &jac_met[qid]);
                let coeff = -factor * quad_w[qid] * tau * jac_det[qid];
                for v in 0..num_node {
                    for j in 0..num_node {
                        let vel_grad_j = vel_x * quad_gnx[qid][j] + vel_y * quad_gny[qid][j];
                        ax_loc[v][j] += coeff * quad_gnx[qid][v] * den * vel_grad_j;
                        ay_loc[v][j] += coeff * quad_gny[qid][v] * den * vel_grad_j;
                        ap_loc[v][j] += coeff * (quad_gnx[qid][v] * quad_gnx[qid][j] + quad_gny[qid][v] * quad_gny[qid][j]);
                    }
                    b_loc[v] += coeff * (quad_gnx[qid][v] * fce_x + quad_gny[qid][v] * fce_y);
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
