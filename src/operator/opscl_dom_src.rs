use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub src_id: usize, // source scalar
    pub unk_id: usize, // unknown scalar
}

impl OpSclDomSource {
    pub fn new(dom_id: usize, src_id: usize, unk_id: usize) -> OpSclDomSource {
        // adds the diffusion term to scalar transport equations
        // d(m_i * c_i)/dt = -div( N_i ) + R_i
        // 
        // src - source scalar (R_i)
        // unk - unknown scalar (c_i)

        // create struct
        let mut oper_src = OpSclDomSource::default();
        oper_src.dom_id = dom_id;
        oper_src.src_id = src_id;
        oper_src.unk_id = unk_id;

        // result
        oper_src
    }
}

impl OperatorBase for OpSclDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the source term
        // +(R, w)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add +(R, w)_dom to A -> add -(R, w)_dom to b
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let src_scl = &vars.scl_dom[self.src_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local vector

            // initialize local vectors
            let num_node = dom.elem_node[eid];
            let mut b_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local vector
            for qid in 0..num_quad {
                let src = src_scl.compute_quad(vars, eid, qid, t);
                let coeff = -factor * quad_w[qid] * src * jac_det[qid];
                for v in 0..num_node {
                    b_loc[v] += coeff * quad_n[qid][v];
                }
            }

            // step 2: add to global matrix

            // iterate over local vector entries
            let node_id = &dom.elem_node_id[eid];
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
