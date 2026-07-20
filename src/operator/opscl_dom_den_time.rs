use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDensityTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub den_id: usize, // density
    pub unk_id: usize, // unknown scalar (pressure)
}

impl OpSclDomDensityTime {
    pub fn new(dom_id: usize, den_id: usize, unk_id: usize) -> OpSclDomDensityTime {
        // adds the density time derivative to the continuity equation
        // d(rho)/dt = -div(rho * v)
        //
        // den - density (rho)
        // unk - unknown scalar (equation added to rows of this scalar; e.g., pressure)

        // create struct
        let mut oper_time = OpSclDomDensityTime::default();
        oper_time.dom_id = dom_id;
        oper_time.den_id = den_id;
        oper_time.unk_id = unk_id;

        // result
        oper_time
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t_next: f64, dt: f64, factor: f64) {
        // density is evaluated from the current nonlinear iterate and the
        // converged solution at the previous time step
        // d(rho)/dt = (rho_next - rho_curr) / dt
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // add +(d(rho)/dt, w)_dom to b

        // get objects
        let dom = &vars.dom[self.dom_id];
        let itgdom = &vars.itg_dom[self.dom_id];
        let den_scl = &vars.scl_dom[self.den_id];
        let unk_scl = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local vector

            // initialize local vector
            let num_node = dom.elem_node[eid];
            let mut b_loc = vec![0.0; num_node];

            // get quadrature point data
            let num_quad = itgdom.num_quad[eid];
            let quad_w = &itgdom.quad_w[eid];
            let quad_n = &itgdom.quad_n[eid];
            let jac_det = &itgdom.jac_det[eid];

            // assemble local vector
            for qid in 0..num_quad {
                let den_next = den_scl.compute_quad(vars, eid, qid, t_next);
                let t_curr = t_next - dt;
                let den_curr = den_scl.compute_quad_prev(vars, eid, qid, t_curr);
                let coeff = factor * quad_w[qid] * (den_next - den_curr) * jac_det[qid] / dt;
                for v in 0..num_node {
                    b_loc[v] += coeff * quad_n[qid][v];
                }
            }

            // step 2: add to global vector
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

impl OperatorBase for OpSclDomDensityTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, _factor: f64) {
        // do not use this. use apply_time instead
        panic!("Used apply for OpSclDomDensityTime. Must use apply_time instead.");
    }
}
