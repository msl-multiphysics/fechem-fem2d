use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomDiffusion {
    // domain
    pub dom_id: usize,

    // scalars
    pub visc_id: usize, // viscosity
    pub unk_id: usize, // unknown vector
    pub drv_id: usize, // driving vector
}

impl OpVecDomDiffusion {
    pub fn new(dom_id: usize, visc_id: usize, unk_id: usize, drv_id: usize) -> OpVecDomDiffusion {
        // adds the diffusion term to the momentum transport equation
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        // T_i += -mu * (grad(v) + grad(v)^T) + (2/3) * mu * div(v) * I
        // 
        // visc - viscosity (mu)
        // unk - unknown vector (v_i)
        // drv - driving vector (v_j)
        
        // create struct
        let mut oper_diff = OpVecDomDiffusion::default();
        oper_diff.dom_id = dom_id;
        oper_diff.visc_id = visc_id;
        oper_diff.drv_id = drv_id;
        oper_diff.unk_id = unk_id;

        // result
        oper_diff
    }
}

impl OperatorBase for OpVecDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the diffusion term
        // +(div(mu * (grad(v) + grad(v)^T) - (2/3) * mu * div(v) * I), w)_dom
        // = -(mu * (grad(v) + grad(v)^T), grad(w))_dom + ((2/3) * mu * div(v), div(w))_dom + (tau . n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(mu * (grad(v) + grad(v)^T), grad(w))_dom + ((2/3) * mu * div(v), div(w))_dom to A
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let visc_scl = &vars.scl_dom[self.visc_id];
        let unk_scl = &vars.vec_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut axx_loc = vec![vec![0.0; num_node]; num_node]; // x momentum, x velocity
            let mut axy_loc = vec![vec![0.0; num_node]; num_node]; // x momentum, y velocity
            let mut ayx_loc = vec![vec![0.0; num_node]; num_node]; // y momentum, x velocity
            let mut ayy_loc = vec![vec![0.0; num_node]; num_node]; // y momentum, y velocity

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let visc = visc_scl.compute_quad(vars, eid, qid, t);
                let coeff = -factor * quad_w[qid] * visc * jac_det[qid];
                let coeff_dil = -(2.0 / 3.0) * coeff; // +(2/3) * mu * (div(v), div(w))
                for v in 0..num_node {
                    for j in 0..num_node {
                        let gnx_v = quad_gnx[qid][v];
                        let gny_v = quad_gny[qid][v];
                        let gnx_j = quad_gnx[qid][j];
                        let gny_j = quad_gny[qid][j];

                        // laplacian: -(mu * grad(v), grad(w))
                        axx_loc[v][j] += coeff * (gnx_v * gnx_j + gny_v * gny_j);
                        ayy_loc[v][j] += coeff * (gnx_v * gnx_j + gny_v * gny_j);

                        // transpose viscosity: -(mu * grad(v)^T, grad(w))
                        axx_loc[v][j] += coeff * gnx_v * gnx_j;
                        axy_loc[v][j] += coeff * gny_v * gnx_j;
                        ayx_loc[v][j] += coeff * gnx_v * gny_j;
                        ayy_loc[v][j] += coeff * gny_v * gny_j;

                        // dilatatory: +((2/3) * mu * div(v), div(w))
                        axx_loc[v][j] += coeff_dil * gnx_v * gnx_j;
                        axy_loc[v][j] += coeff_dil * gnx_v * gny_j;
                        ayx_loc[v][j] += coeff_dil * gny_v * gnx_j;
                        ayy_loc[v][j] += coeff_dil * gny_v * gny_j;
                    }
                }
            }

            // step 2: add to global matrix

            // iterate over local matrix entries
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk_scl.node_dir[nid_v] {
                    continue;
                }
                
                // add to global matrix
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.drv_id, 0, nid_j, axx_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.drv_id, 1, nid_j, axy_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.drv_id, 0, nid_j, ayx_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.drv_id, 1, nid_j, ayy_loc[v][j]);
                }
            }

        }
    }
}
