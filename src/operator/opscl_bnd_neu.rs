use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndNeumann {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub flx_id: usize, // flux
    pub unk_id: usize, // unknown scalar
}

impl OpSclBndNeumann {
    pub fn new(bnd_id: usize, flx_id: usize, unk_id: usize) -> OpSclBndNeumann {
        // applies a Neumann BC
        // d(m_i * c_i)/dt = -div(N_i) + R_i
        // N_i * n = F_i
        // 
        // flx - prescribed normal outward flux (F_i)
        // unk - unknown scalar (c_i)
        // n is the outward normal from the boundary
        // N_i is the flux of c_i (due to diffusion, advection, etc.)

        // create struct
        let mut oper_diff = OpSclBndNeumann::default();
        oper_diff.bnd_id = bnd_id;
        oper_diff.flx_id = flx_id;
        oper_diff.unk_id = unk_id;

        // result
        oper_diff
    }
}

impl OperatorBase for OpSclBndNeumann {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the diffusion term
        // -(div(N, w)_dom = +(N, grad(w))_dom - (N * n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(N * n, w)_bnd to A -> add -(F, w)_bnd to A -> add +(F, w)_bnd to b
        
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let itgbnd = &vars.itg_bnd[self.bnd_id];
        let flx_scl = &vars.scl_bnd[self.flx_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // step 1: assemble local vector

            // initialize local vectors
            let num_node = bnd.elem_node[eid];
            let mut b_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgbnd.num_quad[eid];
            let quad_w = &itgbnd.quad_w[eid];
            let quad_n = &itgbnd.quad_n[eid];
            let jac_det = &itgbnd.jac_det[eid];

            // assemble local vector
            for qid in 0..num_quad {
                let flx = flx_scl.compute_quad(vars, eid, qid, t);
                let coeff = factor * quad_w[qid] * flx * jac_det[qid].sqrt();
                for v in 0..num_node {
                    b_loc[v] += coeff * quad_n[qid][v];
                }
            }

            // step 2: add to global matrix

            // iterate over nodes in element
            let node_id = &bnd.elem_node_id[eid];
            for v in 0..num_node {
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk_scl.node_dir[nid_v] {
                    continue;
                }
                
                // add to global vector
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, b_loc[v]);
            }

        }
    }
}
