use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndPressure {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub pres_id: usize, // prescribed pressure
    pub unk_id: usize,  // unknown vector
}

impl OpVecBndPressure {
    pub fn new(bnd_id: usize, pres_id: usize, unk_id: usize) -> OpVecBndPressure {
        // applies a pressure BC to the momentum transport equation
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        // T_i += -grad(p)
        //
        // pres - prescribed pressure (p)
        // unk - unknown vector (v_i)
        // n is the outward normal from the boundary

        // create struct
        let mut oper_pres = OpVecBndPressure::default();
        oper_pres.bnd_id = bnd_id;
        oper_pres.pres_id = pres_id;
        oper_pres.unk_id = unk_id;

        // result
        oper_pres
    }
}

impl OperatorBase for OpVecBndPressure {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the boundary term from the weak form of the pressure gradient
        // -(grad(p), w)_dom = +(p, div(w))_dom - (pw . n)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(pw . n)_bnd to A -> add +(pw . n)_bnd to b

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let itgbnd = &vars.itg_bnd[self.bnd_id];
        let pres_scl = &vars.scl_bnd[self.pres_id];
        let unk_vec = &vars.vec_dom[self.unk_id];

        // iterate over elements
        for eid in 0..bnd.num_elem {
            // step 1: assemble local vectors

            // initialize local vectors
            let num_node = bnd.elem_node[eid];
            let mut bx_loc = vec![0.0; num_node]; // x momentum
            let mut by_loc = vec![0.0; num_node]; // y momentum

            // get quadrature point data
            let num_quad = itgbnd.num_quad[eid];
            let quad_w = &itgbnd.quad_w[eid];
            let quad_n = &itgbnd.quad_n[eid];
            let jac_mat = &itgbnd.jac_mat[eid];

            // assemble local vectors
            for qid in 0..num_quad {
                let pres = pres_scl.compute_quad(vars, eid, qid, t);
                let coeff = factor * quad_w[qid] * pres;
                // reverse the inward normal stored in the jacobian
                let normal_x = -jac_mat[qid][0][1];
                let normal_y = -jac_mat[qid][1][1];
                for v in 0..num_node {
                    bx_loc[v] += coeff * normal_x * quad_n[qid][v];
                    by_loc[v] += coeff * normal_y * quad_n[qid][v];
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
                if unk_vec.node_dir[nid_dom] {
                    continue;
                }

                // add to global vector
                self.add_b_vecdom(vars, b_vec, self.unk_id, 0, nid_dom, bx_loc[v]);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 1, nid_dom, by_loc[v]);
            }
        }
    }
}
