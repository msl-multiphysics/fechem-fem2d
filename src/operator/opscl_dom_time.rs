use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub wgt_id: usize, // mass scalar
    pub unk_id: usize, // unknown scalar
}

impl OpSclDomTime {
    pub fn new(dom_id: usize, wgt_id: usize, unk_id: usize) -> OpSclDomTime {
        // adds the time derivative term to scalar transport equations
        // d(m_i * c_i)/dt = -div(N_i) + R_i
        // 
        // wgt - mass scalar (m_i)
        // unk - unknown scalar (c_i)

        // create struct
        let mut oper_time = OpSclDomTime::default();
        oper_time.dom_id = dom_id;
        oper_time.wgt_id = wgt_id;
        oper_time.unk_id = unk_id;

        // result
        oper_time
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t_next: f64, dt: f64, factor: f64) {
        // time derivative is discretized using backward Euler
        // d(m * c)/dt = (m_next * c_next - m_curr * c_curr) / dt
        // 
        // apply the weak form of the time derivative term
        // (d(m * c)/dt, w)_dom = ((m_next * c_next)/dt, w)_dom - ((m_curr * c_curr)/dt, w)_dom
        // 
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add ((m_next * c_next)/dt, w)_dom to b -> add -((m_next * c_next)/dt, w)_dom to b
        // add -((m_curr * c_curr)/dt, w)_dom to b
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let wgt_scl = &vars.scl_dom[self.wgt_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

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
                let wgt = wgt_scl.compute_quad(vars, eid, qid, t_next);
                let coeff = -factor * quad_w[qid] * wgt * jac_det[qid] / dt;

                // current time step
                let t_curr = t_next - dt;
                let wgt_curr = wgt_scl.compute_quad_prev(vars, eid, qid, t_curr);
                let coeff_curr = -factor * quad_w[qid] * wgt_curr * jac_det[qid] / dt;

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
                if unk_scl.node_dir[nid_v] {
                    continue;
                }

                // add next time step
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    self.add_a_scldom(vars, a_triplet, self.unk_id, nid_v, self.unk_id, nid_j, an_loc[v][j]);
                }

                // add current time step
                let mut ac_sum = 0.0;
                for j in 0..num_node {
                    let nid_j = node_id[j];
                    ac_sum += ac_loc[v][j] * unk_scl.node_value_prev[nid_j];
                }
                self.add_b_scldom(vars, b_vec, self.unk_id, nid_v, ac_sum);
            }

        }
    }
}

impl OperatorBase for OpSclDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, _factor: f64) {
        // do not use this. use apply_time instead
        panic!("Used apply for OperatorTime. Must use apply_time instead.");
    }
}
