use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndDivergence {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub den_id: usize, // density
    pub vel_id: usize, // prescribed velocity
    pub unk_id: usize, // unknown scalar
}

impl OpSclBndDivergence {
    pub fn new(bnd_id: usize, den_id: usize, vel_id: usize, unk_id: usize) -> OpSclBndDivergence {
        // applies a velocity BC to the continuity equation
        // d(rho)/dt = -div(rho * v)
        //
        // den - density (rho)
        // vel - prescribed velocity (v)
        // unk - unknown scalar (equation added to rows of this scalar; e.g., pressure)
        // n is the outward normal from the boundary

        // create struct
        let mut oper_div = OpSclBndDivergence::default();
        oper_div.bnd_id = bnd_id;
        oper_div.den_id = den_id;
        oper_div.vel_id = vel_id;
        oper_div.unk_id = unk_id;

        // result
        oper_div
    }
}

impl OperatorBase for OpSclBndDivergence {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the boundary term from the weak form of the divergence term
        // -(div(rho * v), w)_dom = +(rho * v, grad(w))_dom - (rho * v . n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(rho * v . n, w)_bnd to A -> add +(rho * v . n, w)_bnd to b

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let itgbnd = &vars.itg_bnd[self.bnd_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let vel_vec = &vars.vec_bnd[self.vel_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // step 1: assemble local vector

            // initialize local vector
            let num_node = bnd.elem_node[eid];
            let mut b_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgbnd.num_quad[eid];
            let quad_w = &itgbnd.quad_w[eid];
            let quad_n = &itgbnd.quad_n[eid];
            let jac_mat = &itgbnd.jac_mat[eid];

            // assemble local vector
            for qid in 0..num_quad {
                let den = den_scl.compute_quad_unknown_boundary(bnd, itgbnd, eid, qid);
                let (vel_x, vel_y) = vel_vec.compute_quad(vars, eid, qid, t);
                // reverse the inward normal stored in the jacobian
                let normal_x = -jac_mat[qid][0][1];
                let normal_y = -jac_mat[qid][1][1];
                let flux = den * (vel_x * normal_x + vel_y * normal_y);
                let coeff = factor * quad_w[qid] * flux;
                for v in 0..num_node {
                    b_loc[v] += coeff * quad_n[qid][v];
                }
            }

            // step 2: add to global vector

            // iterate over nodes in element
            let node_id = &bnd.elem_node_id[eid];
            for v in 0..num_node {
                // get node ids
                let nid_bnd = node_id[v];
                let nid_dom = bnd.node_bnd_dom_id[nid_bnd];

                // skip if dirichlet BC
                if unk_scl.node_dir[nid_dom] {
                    continue;
                }

                // add to global vector
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_dom, b_loc[v]);
            }
        }
    }
}
