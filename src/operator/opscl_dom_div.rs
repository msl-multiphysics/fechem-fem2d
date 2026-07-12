use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDivergence {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub vel_id: usize, // velocity
    pub unk_id: usize, // unknown scalar
}

impl OpSclDomDivergence {
    pub fn new(dom_id: usize, den_id: usize, vel_id: usize, unk_id: usize) -> OpSclDomDivergence {
        // adds the divergence term to the continuity equation
        // d(rho)/dt = -div(rho * v)
        // 
        // den - density (rho)
        // vel - velocity (v)
        // unk - unknown scalar (equation added to rows of this scalar; e.g., pressure)
        
        // create struct
        let mut oper_diff = OpSclDomDivergence::default();
        oper_diff.dom_id = dom_id;
        oper_diff.den_id = den_id;
        oper_diff.vel_id = vel_id;
        oper_diff.unk_id = unk_id;

        // result
        oper_diff
    }
}

impl OperatorBase for OpSclDomDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the divergence term
        // -(div(rho * v), w)_dom = +(rho * v, grad(w))_dom - (rho * v * n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add +(rho * v, grad(w))_dom to A
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut ax_loc = vec![vec![0.0; num_node]; num_node];  // x-component of velocity
            let mut ay_loc = vec![vec![0.0; num_node]; num_node];  // y-component of velocity

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let den = den_scl.compute_quad(vars, eid, qid, t);
                let coeff = factor * quad_w[qid] * den * jac_det[qid];
                for v in 0..num_node {
                    for j in 0..num_node {
                        ax_loc[v][j] += coeff * quad_n[qid][j] * quad_gnx[qid][v];
                        ay_loc[v][j] += coeff * quad_n[qid][j] * quad_gny[qid][v];
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
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 0, nid_j, ax_loc[v][j]);
                    self.add_a_scldom_vecdom(vars, a_triplet, self.unk_id, nid_v, self.vel_id, 1, nid_j, ay_loc[v][j]);
                }
            }

        }
    }
}
