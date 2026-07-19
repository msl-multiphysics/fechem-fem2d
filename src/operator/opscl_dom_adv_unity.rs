use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomAdvectionUnity {
    // domain
    pub dom_id: usize,

    // scalars
    pub vel_id: usize, // velocity vector
    pub unk_id: usize, // unknown scalar
}

impl OpSclDomAdvectionUnity {
    pub fn new(dom_id: usize, vel_id: usize, unk_id: usize) -> OpSclDomAdvectionUnity {
        // adds the advective term to scalar transport equations
        // dc_i/dt = -div(N_i) + R_i
        // N_i += c_i * v
        // 
        // vel - velocity vector (v)
        // unk - unknown scalar (c_i)
        // weight is unity (1)

        // create struct
        let mut oper_adv = OpSclDomAdvectionUnity::default();
        oper_adv.dom_id = dom_id;
        oper_adv.vel_id = vel_id;
        oper_adv.unk_id = unk_id;

        // result
        oper_adv
    }
}

impl OperatorBase for OpSclDomAdvectionUnity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the advective term
        // -(div(c * v), w)_dom = +(c * v, grad(w))_dom - (c * v . n, w)_bnd
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add +(c * v, grad(w))_dom to A
        // the boundary term is assembled separately by OpSclBndOutflow
        // v is lagged by 1 iteration
    
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let vel_vec = &vars.vec_dom[self.vel_id];
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
            let quad_n = &itgdom.quad_n[eid];
            let quad_gnx = &itgdom.quad_gnx[eid];
            let quad_gny = &itgdom.quad_gny[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let (vel_x, vel_y) = vel_vec.compute_quad(vars, eid, qid, t);  // lag the velocity by 1 iteration
                let coeff = factor * quad_w[qid] * jac_det[qid];  // weight = 1
                for v in 0..num_node {
                    for j in 0..num_node {
                        a_loc[v][j] += coeff * (vel_x * quad_gnx[qid][v] + vel_y * quad_gny[qid][v]) * quad_n[qid][j];
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
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, self.unk_id, nid_j, a_loc[v][j]);
                }
            }

        }
    }
}
