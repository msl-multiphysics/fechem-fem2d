use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub src_id: usize, // source
    pub unk_id: usize, // unknown vector

}

impl OpVecDomSource {
    pub fn new(dom_id: usize, src_id: usize, unk_id: usize) -> OpVecDomSource {
        // adds the source term to the momentum transport equation
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        // 
        // src - source vector (f_i)
        // unk - unknown vector (v_i)

        // create struct
        let mut oper_src = OpVecDomSource::default();
        oper_src.dom_id = dom_id;
        oper_src.src_id = src_id;
        oper_src.unk_id = unk_id;

        // result
        oper_src
    }
}

impl OperatorBase for OpVecDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the source term
        // +(f, w)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add +(f, w)_dom to A -> add -(f, w)_dom to b
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let src_vec = &vars.vec_dom[self.src_id];
        let unk_vec = &vars.vec_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut bx_loc = vec![0.0; num_node];
            let mut by_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let (src_x, src_y) = src_vec.compute_quad(vars, eid, qid, t);
                let coeff = -factor * quad_w[qid] * jac_det[qid];
                for v in 0..num_node {
                    bx_loc[v] += coeff * src_x * quad_n[qid][v];
                    by_loc[v] += coeff * src_y * quad_n[qid][v];
                }
            }

            // step 2: assemble global matrix

            // iterate over local matrix entries
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                // skip if dirichlet BC
                let nid_v = node_id[v];
                if unk_vec.node_dir[nid_v] {
                    continue;
                }
                
                // add to global matrix
                self.add_b_vecdom(vars, b_vec, self.unk_id, 0, nid_v, bx_loc[v]);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 1, nid_v, by_loc[v]);
            }
        }
    }
}
