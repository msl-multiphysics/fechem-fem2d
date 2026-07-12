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
        // T_i += -mu_ij * grad(v_j)
        // 
        // visc - viscosity (mu_ij)
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
        // +(div(mu * grad(v), w)_dom = -(mu * grad(v), grad(w))_dom + (mu * grad(v) * n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(mu * grad(v), grad(w))_dom to A
    
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
            let mut a_loc = vec![vec![0.0; num_node]; num_node];  // both x and y momentum have the same local matrix

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
                for v in 0..num_node {
                    for j in 0..num_node {
                        a_loc[v][j] += coeff * (quad_gnx[qid][v] * quad_gnx[qid][j] + quad_gny[qid][v] * quad_gny[qid][j]);
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
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.drv_id, 0, nid_j, a_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.drv_id, 1, nid_j, a_loc[v][j]);
                }
            }

        }
    }
}
