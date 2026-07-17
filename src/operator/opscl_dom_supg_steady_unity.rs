use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomSupgSteadyUnity {
    // domain
    pub dom_id: usize,

    // scalars
    pub tau_diff_id: usize, // effective diffusion coefficient used in tau
    pub src_id: usize, // source scalar
    pub unk_id: usize, // unknown scalar
    pub diff_drv_ids: Vec<(usize, usize)>, // (diffusion coefficient, driving scalar)

    // vectors
    pub vel_id: usize, // velocity vector
}

impl OpSclDomSupgSteadyUnity {
    pub fn new(dom_id: usize, tau_diff_id: usize, vel_id: usize, src_id: usize, unk_id: usize, diff_drv_ids: Vec<(usize, usize)>) -> OpSclDomSupgSteadyUnity {
        // adds steady SUPG stabilization to a scalar transport equation
        // dc_i/dt = -div(c_i * v - sum_j(D_ij * grad(c_j))) + R_i
        //
        // tau_diff - effective diffusion coefficient used in tau (normally D_ii)
        // vel - velocity vector (v)
        // src - source scalar (R_i)
        // unk - unknown scalar (c_i)
        // diff_drv_ids - (D_ij, c_j) pairs in the diffusive flux
        // weight is unity (1)

        // create struct
        let mut oper_supg = OpSclDomSupgSteadyUnity::default();
        oper_supg.dom_id = dom_id;
        oper_supg.tau_diff_id = tau_diff_id;
        oper_supg.vel_id = vel_id;
        oper_supg.src_id = src_id;
        oper_supg.unk_id = unk_id;
        oper_supg.diff_drv_ids = diff_drv_ids;

        // result
        oper_supg
    }

    fn compute_tau(&self, diff_val: f64, vel_x: f64, vel_y: f64, jac_met: &[[f64; 2]; 2]) -> f64 {
        // metric-based steady SUPG time scale
        let g00 = jac_met[0][0];
        let g01 = jac_met[0][1];
        let g10 = jac_met[1][0];
        let g11 = jac_met[1][1];
        let adv = (vel_x * (g00 * vel_x + g01 * vel_y) + vel_y * (g10 * vel_x + g11 * vel_y)).max(0.0).sqrt();
        let diff = diff_val.abs() * (g00 * g00 + g01 * g01 + g10 * g10 + g11 * g11).sqrt();

        1.0 / (2.0 * adv + 4.0 * diff + 1e-30)
    }
}

impl OperatorBase for OpSclDomSupgSteadyUnity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the steady SUPG stabilization term
        // tau * (v . grad(w), div(c * v) - div(sum_j(D_ij * grad(c_j))) - R_i)_dom
        // on P1 elements, div(D_ij * grad(c_j)) reduces to grad(D_ij) . grad(c_j)
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add the SUPG stabilization contributions to A and b

        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let tau_diff_scl = &vars.scl_dom[self.tau_diff_id];
        let src_scl = &vars.scl_dom[self.src_id];
        let unk_scl = &vars.scl_dom[self.unk_id];
        let vel_vec = &vars.vec_dom[self.vel_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrices and vector

            // initialize local matrices and vector
            let num_node = dom.elem_node[eid];
            let mut adv_loc = vec![vec![0.0; num_node]; num_node];
            let mut diff_loc = vec![vec![vec![0.0; num_node]; num_node]; self.diff_drv_ids.len()];
            let mut b_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];
            let jac_met = &itgdom.jac_met[eid];

            // assemble local matrices and vector
            for qid in 0..num_quad {
                let tau_diff = tau_diff_scl.compute_quad(vars, eid, qid, t);
                let src = src_scl.compute_quad(vars, eid, qid, t);
                let (vel_x, vel_y) = vel_vec.compute_quad(vars, eid, qid, t);  // lag the velocity by 1 iteration
                let vel_grad = vel_vec.compute_quad_grad(vars, eid, qid, t);
                let div_vel = vel_grad[0][0] + vel_grad[1][1];
                let tau = self.compute_tau(tau_diff, vel_x, vel_y, &jac_met[qid]);
                let coeff = -factor * quad_w[qid] * tau * jac_det[qid];

                for v in 0..num_node {
                    let vel_grad_v = vel_x * quad_gnx[qid][v] + vel_y * quad_gny[qid][v];
                    for j in 0..num_node {
                        let vel_grad_j = vel_x * quad_gnx[qid][j] + vel_y * quad_gny[qid][j];
                        let div_adv_j = vel_grad_j + div_vel * quad_n[qid][j];
                        adv_loc[v][j] += coeff * vel_grad_v * div_adv_j;
                    }
                    b_loc[v] += coeff * vel_grad_v * src;
                }

                for (did, &(diff_id, _)) in self.diff_drv_ids.iter().enumerate() {
                    let [diff_x, diff_y] = vars.scl_dom[diff_id].compute_quad_grad(vars, eid, qid, t);
                    for v in 0..num_node {
                        let vel_grad_v = vel_x * quad_gnx[qid][v] + vel_y * quad_gny[qid][v];
                        for j in 0..num_node {
                            let div_diff_j = diff_x * quad_gnx[qid][j] + diff_y * quad_gny[qid][j];
                            diff_loc[did][v][j] += -coeff * vel_grad_v * div_diff_j;
                        }
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

                // add advection and diffusion matrices
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, self.unk_id, nid_j, adv_loc[v][j]);
                    for (did, &(_, drv_id)) in self.diff_drv_ids.iter().enumerate() {
                        self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, drv_id, nid_j, diff_loc[did][v][j]);
                    }
                }

                // add source vector
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, b_loc[v]);
            }
        }
    }
}
