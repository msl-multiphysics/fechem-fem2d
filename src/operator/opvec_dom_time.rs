use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub unk_id: usize, // unknown vector
}

impl OpVecDomTime {
    pub fn new(dom_id: usize, den_id: usize, unk_id: usize) -> OpVecDomTime {
        // adds the time derivative term to vector transport equations
        // d(den_i * v_i)/dt = -div(T_i) + f_i
        // 
        // den - density (den_i)
        // unk - unknown vector (v_i)

        // create struct
        let mut oper_time = OpVecDomTime::default();
        oper_time.dom_id = dom_id;
        oper_time.den_id = den_id;
        oper_time.unk_id = unk_id;

        // result
        oper_time
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t_next: f64, dt: f64, factor: f64) {
        // time derivative is discretized using backward Euler
        // d(den * v)/dt = (den_next * v_next - den_curr * v_curr) / dt
        // 
        // apply the weak form of the time derivative term
        // (d(den * v)/dt, w)_dom = ((den_next * v_next)/dt, w)_dom - ((den_curr * v_curr)/dt, w)_dom
        // 
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add ((den_next * v_next)/dt, w)_dom to b -> add -((den_next * v_next)/dt, w)_dom to b
        // add -((den_curr * v_curr)/dt, w)_dom to b
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let unk_vec = &vars.vec_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local mass matrix

            // initialize local matrix
            let num_node = dom.elem_node[eid];
            let mut an_loc = vec![vec![0.0; num_node]; num_node];  // next time step
            let mut ac_loc = vec![vec![0.0; num_node]; num_node];  // current time step

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local mass matrix
            for qid in 0..num_quad {
                // next time step
                let den = den_scl.compute_quad(vars, eid, qid, t_next);
                let coeff = -factor * quad_w[qid] * den * jac_det[qid] / dt;

                // current time step
                let t_curr = t_next - dt;
                let den_curr = den_scl.compute_quad_prev(vars, eid, qid, t_curr);
                let coeff_curr = -factor * quad_w[qid] * den_curr * jac_det[qid] / dt;

                // load entries
                for v in 0..num_node {
                    for j in 0..num_node {
                        an_loc[v][j] += coeff * quad_n[qid][v] * quad_n[qid][j];
                        ac_loc[v][j] += coeff_curr * quad_n[qid][v] * quad_n[qid][j];
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

                // add next time step
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 0, nid_v, self.unk_id, 0, nid_j, an_loc[v][j]);
                    self.add_a_vecdom(vars, a_triplet, self.unk_id, 1, nid_v, self.unk_id, 1, nid_j, an_loc[v][j]);
                }

                // add current time step
                let mut ac_sum_x = 0.0;
                let mut ac_sum_y = 0.0;
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    ac_sum_x += ac_loc[v][j] * unk_vec.node_value_prev_x[nid_j];
                    ac_sum_y += ac_loc[v][j] * unk_vec.node_value_prev_y[nid_j];
                }
                self.add_b_vecdom(vars, b_vec, self.unk_id, 0, nid_v, ac_sum_x);
                self.add_b_vecdom(vars, b_vec, self.unk_id, 1, nid_v, ac_sum_y);
            }

        }
    }
}

impl OperatorBase for OpVecDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, _factor: f64) {
        // do not use this. use apply_time instead
        panic!("Used apply for OperatorTime. Must use apply_time instead.");
    }
}
