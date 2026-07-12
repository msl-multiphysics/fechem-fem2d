use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
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
        // adds SUPG stabilization to the momentum transport equation
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        // 
        // den - density (den_i)
        // visc - viscosity (mu)
        // unk - unknown vector (v_i)
        // pres - pressure (p)
        // fce - body force (f)
        // drv - driving vector (v_j)

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
        // applies the weak form of the SUPG stabilization term
        // tau * (v_j . grad(w), rho * v_j . grad(v_i) + grad_i(p) - f_i)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add the SUPG stabilization contributions to A and b
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let visc_scl = &vars.scl_dom[self.visc_id];
        let fce_vec = &vars.vec_dom[self.fce_id];
        let unk_vec = &vars.vec_dom[self.unk_id];
        let drv_vec = &vars.vec_dom[self.drv_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix and vector

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];  // both x and y momentum have the same local matrix
            let mut axp_loc = vec![vec![0.0; num_node]; num_node];  // pressure coupling in x momentum
            let mut ayp_loc = vec![vec![0.0; num_node]; num_node];  // pressure coupling in y momentum
            let mut bx_loc = vec![0.0; num_node];
            let mut by_loc = vec![0.0; num_node];

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
                let (drv_x, drv_y) = drv_vec.compute_quad(vars, eid, qid, t);  // lag the driving vector by 1 iteration
                let (fce_x, fce_y) = fce_vec.compute_quad(vars, eid, qid, t);
                let tau = self.compute_tau(den, visc, drv_x, drv_y, &jac_met[qid]);
                let coeff = -factor * quad_w[qid] * tau * jac_det[qid];
                for v in 0..num_node {
                    let drv_grad_v = drv_x * quad_gnx[qid][v] + drv_y * quad_gny[qid][v];
                    for j in 0..num_node {
                        let drv_grad_j = drv_x * quad_gnx[qid][j] + drv_y * quad_gny[qid][j];
                        a_loc[v][j] += coeff * drv_grad_v * den * drv_grad_j;
                        axp_loc[v][j] += coeff * drv_grad_v * quad_gnx[qid][j];
                        ayp_loc[v][j] += coeff * drv_grad_v * quad_gny[qid][j];
                    }
                    bx_loc[v] += coeff * drv_grad_v * fce_x;
                    by_loc[v] += coeff * drv_grad_v * fce_y;
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
