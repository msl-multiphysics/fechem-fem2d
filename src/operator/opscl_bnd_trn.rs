use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndTransfer {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub trn_id: usize, // transfer coefficient
    pub ext_id: usize, // external value
    pub unk_id: usize, // unknown scalar
}

impl OpSclBndTransfer {
    pub fn new(bnd_id: usize, trn_id: usize, ext_id: usize, unk_id: usize) -> OpSclBndTransfer {
        // applies a transfer BC
        // d(m_i * c_i)/dt = -div(N_i) + R_i
        // N_i . n = h_i * (c_i - c_ext)
        // 
        // trn - transfer coefficient (h_i)
        // ext - external value (c_ext)
        // unk - unknown scalar (c_i)
        // n is the outward normal from the boundary
        // N_i is the flux of c_i (due to diffusion, advection, etc.)

        // create struct
        let mut oper_trn = OpSclBndTransfer::default();
        oper_trn.bnd_id = bnd_id;
        oper_trn.trn_id = trn_id;
        oper_trn.ext_id = ext_id;
        oper_trn.unk_id = unk_id;

        // result
        oper_trn
    }
}

impl OperatorBase for OpSclBndTransfer {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the transfer term
        // -(div(N), w)_dom = +(N, grad(w))_dom - (N . n, w)_bnd
        // with N . n = h * (c - c_ext)
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(N . n, w)_bnd to A -> add -(h * (c - c_ext), w)_bnd to A
        // add -(h * c, w)_bnd to A and -(h * c_ext, w)_bnd to b
        
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let itgbnd = &vars.itg_bnd[self.bnd_id];
        let trn_scl = &vars.scl_bnd[self.trn_id];
        let ext_scl = &vars.scl_bnd[self.ext_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // step 1: assemble local matrix and vector

            // initialize local matrix and vector
            let num_node = bnd.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];
            let mut b_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgbnd.num_quad[eid];
            let quad_w = &itgbnd.quad_w[eid];
            let quad_n = &itgbnd.quad_n[eid];
            let jac_det = &itgbnd.jac_det[eid];

            // assemble local matrix and vector
            for qid in 0..num_quad {
                let trn = trn_scl.compute_quad(vars, eid, qid, t);
                let ext = ext_scl.compute_quad(vars, eid, qid, t);
                let coeff = -factor * quad_w[qid] * jac_det[qid].sqrt();
                for v in 0..num_node {
                    for j in 0..num_node {
                        a_loc[v][j] += coeff * trn * quad_n[qid][v] * quad_n[qid][j];
                    }
                    b_loc[v] += coeff * trn * ext * quad_n[qid][v];
                }
            }

            // step 2: add to global matrix and vector

            // iterate over nodes in element
            let node_id = &bnd.elem_node_id[eid];
            for v in 0..num_node {
                // get node ids
                let nid_bnd_v = node_id[v];
                let nid_dom_v = bnd.node_bnd_dom_id[nid_bnd_v];

                // skip if dirichlet BC
                if unk_scl.node_dir[nid_dom_v] {
                    continue;
                }
                
                // add to global vector
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_dom_v, b_loc[v]);

                // add to global matrix
                for j in 0..num_node {
                    let nid_bnd_j = node_id[j];
                    let nid_dom_j = bnd.node_bnd_dom_id[nid_bnd_j];
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_dom_v, self.unk_id, nid_dom_j, a_loc[v][j]);
                }
            }

        }
    }
}
