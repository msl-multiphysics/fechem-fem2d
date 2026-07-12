use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomPressure {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize, // unknown vector
    pub pres_id: usize, // pressure
}

impl OpVecDomPressure {
    pub fn new(dom_id: usize, unk_id: usize, pres_id: usize) -> OpVecDomPressure {
        // adds the pressure gradient term to the momentum transport equation
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        // T_i += -grad(p)
        // 
        // unk - unknown vector (v_i)
        // pres - pressure (p)

        // create struct
        let mut oper_pres = OpVecDomPressure::default();
        oper_pres.dom_id = dom_id;
        oper_pres.pres_id = pres_id;
        oper_pres.unk_id = unk_id;

        // result
        oper_pres
    }
}

impl OperatorBase for OpVecDomPressure {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // applies the weak form of the pressure gradient term
        // -(grad(p), w)_dom = +(p, div(w))_dom - (pw . n)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add +(p, div(w))_dom to A
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let unk_vec = &vars.vec_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut ax_loc = vec![vec![0.0; num_node]; num_node];  // x momentum
            let mut ay_loc = vec![vec![0.0; num_node]; num_node];  // y momentum

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let coeff = factor * quad_w[qid] * jac_det[qid];
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
                if unk_vec.node_dir[nid_v] {
                    continue;
                }
                
                // add to global matrix
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_vecdom_scldom(vars, a_triplet, self.unk_id, 0, nid_v, self.pres_id, nid_j, ax_loc[v][j]);
                    self.add_a_vecdom_scldom(vars, a_triplet, self.unk_id, 1, nid_v, self.pres_id, nid_j, ay_loc[v][j]);
                }
            }
        }
    }
}
