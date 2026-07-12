use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomAdvection {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub unk_id: usize, // unknown vector
    pub drv_id: usize, // driving vector
}

impl OpVecDomAdvection {
    pub fn new(dom_id: usize, den_id: usize, unk_id: usize, drv_id: usize) -> OpVecDomAdvection {
        // adds the advective term to the momentum transport equation
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        // T_i += den_i * v_i * v_j
        // 
        // den - density (den_i)
        // unk - unknown vector (v_i)
        // drv - driving vector (v_j)

        // create struct
        let mut oper_adv = OpVecDomAdvection::default();
        oper_adv.dom_id = dom_id;
        oper_adv.den_id = den_id;
        oper_adv.unk_id = unk_id;
        oper_adv.drv_id = drv_id;
        
        // result
        oper_adv
    }
}

impl OperatorBase for OpVecDomAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the advective term
        // -(den * v_j . grad(v_i), w)_dom
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add -(den * v_j . grad(v_i), w)_dom to A
        // v_j is lagged by 1 iteration
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let unk_vec = &vars.vec_dom[self.unk_id];
        let drv_vec = &vars.vec_dom[self.drv_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];  // x and y matrices are the same

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
                let (drv_x, drv_y) = drv_vec.compute_quad(vars, eid, qid, t);  // lag the driving vector by 1 iteration
                let coeff = -factor * quad_w[qid] * den * jac_det[qid];
                for v in 0..num_node {
                    for j in 0..num_node {
                        a_loc[v][j] += coeff * (drv_x * quad_gnx[qid][j] + drv_y * quad_gny[qid][j]) * quad_n[qid][v];
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
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.unk_id, 0, nid_j, a_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.unk_id, 1, nid_j, a_loc[v][j]);
                }
            }
        }
    }
}
