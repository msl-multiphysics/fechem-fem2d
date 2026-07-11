use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDiffusion {
    // domain
    pub dom_id: usize,

    // scalars
    pub diff_id: usize, // diffusion coefficient
    pub unk_id: usize,  // unknown scalar
    pub drv_id: usize,  // driving scalar
}

impl OpSclDomDiffusion {
    pub fn new(dom_id: usize, diff_id: usize, unk_id: usize, drv_id: usize) -> OpSclDomDiffusion {
        // adds the diffusion term to scalar transport equations
        // d(m_i * c_i)/dt = -div( -D_ij * grad(c_j) + v * c_i ) + R_i
        // 
        // diff - diffusion coefficient (D_ij)
        // drv - driving scalar (c_j)
        // unk - unknown scalar (c_i)

        // create struct
        let mut oper_diff = OpSclDomDiffusion::default();
        oper_diff.dom_id = dom_id;
        oper_diff.diff_id = diff_id;
        oper_diff.unk_id = unk_id;
        oper_diff.drv_id = drv_id;

        // result
        oper_diff
    }
}

impl OperatorBase for OpSclDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the diffusion term
        // +(div(D * grad(c), w)_dom = -(D * grad(c), grad(w))_dom + (D * grad(c) * n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(D * grad(c), grad(w))_dom to A
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let diff_scl = &vars.scl_dom[self.diff_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrix
            let num_node = dom.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let diff = diff_scl.compute_quad(vars, eid, qid, t);
                let coeff = -factor * quad_w[qid] * diff * jac_det[qid];
                for v in 0..num_node {
                    for j in 0..num_node {
                        a_loc[v][j] += coeff * (quad_gnx[qid][v] * quad_gnx[qid][j] + quad_gny[qid][v] * quad_gny[qid][j]);
                    }
                }
            }

            // step 2: add to global matrix

            // iterate over nodes in element
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
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, self.drv_id, nid_j, a_loc[v][j]);
                }
            }

        }
    }
}
