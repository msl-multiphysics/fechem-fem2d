use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndOutflow {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub wgt_id: usize, // weighting scalar
    pub vel_id: usize, // velocity vector
    pub unk_id: usize, // unknown scalar
}

impl OpSclBndOutflow {
    pub fn new(bnd_id: usize, wgt_id: usize, vel_id: usize, unk_id: usize) -> OpSclBndOutflow {
        // applies an outflow BC to scalar transport equations
        // d(m_i * c_i)/dt = -div(N_i) + R_i
        // N_i += m_i * c_i * v
        // N_i . n = m_i * c_i * max(v . n, 0)   (outgoing only; zero conductive flux)
        //
        // wgt - weighting scalar (m_i)
        // vel - velocity vector (v)
        // unk - unknown scalar (c_i)
        // n is the outward normal from the boundary
        // completes the boundary term omitted by OpSclDomAdvection
        // max(v . n, 0) suppresses backflow, which would otherwise act as an
        // inflow without an exterior value and create outlet artifacts

        // create struct
        let mut oper_out = OpSclBndOutflow::default();
        oper_out.bnd_id = bnd_id;
        oper_out.wgt_id = wgt_id;
        oper_out.vel_id = vel_id;
        oper_out.unk_id = unk_id;

        // result
        oper_out
    }
}

impl OperatorBase for OpSclBndOutflow {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // applies the boundary term from the weak form of the advective term
        // -(div(m * c * v), w)_dom = +(m * c * v, grad(w))_dom - (m * c * v . n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(m * c * max(v . n, 0), w)_bnd to A
        // v is lagged by 1 iteration

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let itgbnd = &vars.itg_bnd[self.bnd_id];
        let wgt_scl = &vars.scl_dom[self.wgt_id];
        let vel_vec = &vars.vec_dom[self.vel_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // step 1: assemble local matrix

            // initialize local matrix
            let num_node = bnd.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];

            // get quadrature point data
            let num_quad = itgbnd.num_quad[eid];
            let quad_w = &itgbnd.quad_w[eid];
            let quad_n = &itgbnd.quad_n[eid];
            let jac_mat = &itgbnd.jac_mat[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let wgt = wgt_scl.compute_quad_unknown_boundary(bnd, itgbnd, eid, qid);
                let (vel_x, vel_y) = vel_vec.compute_quad_unknown_boundary(bnd, itgbnd, eid, qid);  // lag the velocity by 1 iteration
                // reverse the inward normal stored in the jacobian
                let normal_x = -jac_mat[qid][0][1];
                let normal_y = -jac_mat[qid][1][1];
                let vn = vel_x * normal_x + vel_y * normal_y;
                let vn_out = vn.max(0.0);
                let coeff = -factor * quad_w[qid] * wgt * vn_out;
                for v in 0..num_node {
                    for j in 0..num_node {
                        a_loc[v][j] += coeff * quad_n[qid][v] * quad_n[qid][j];
                    }
                }
            }

            // step 2: add to global matrix

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
