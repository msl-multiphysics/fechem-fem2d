use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OperatorDiffusion {
    // domain
    pub dom_id: usize,

    // scalars
    pub diff_id: usize, // diffusion coefficient
    pub unk_id: usize,  // unknown scalar
    pub drv_id: usize,  // driving scalar
}

impl OperatorDiffusion {
    pub fn new(dom_id: usize, diff_id: usize, unk_id: usize, drv_id: usize) -> OperatorDiffusion {
        // applies -diff * lapl(drv) to the unknown scalar

        // create struct
        let mut oper_diff = OperatorDiffusion::default();
        oper_diff.dom_id = dom_id;
        oper_diff.diff_id = diff_id;
        oper_diff.unk_id = unk_id;
        oper_diff.drv_id = drv_id;

        // result
        oper_diff
    }
}

impl OperatorBase for OperatorDiffusion {
    fn apply(
        &self,
        vars: &Variables,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        _b_vec: &mut Col<f64>,
        t: f64,
        factor: f64,
    ) {
        // get objects
        let dom = &vars.dom[self.dom_id];
        let itg = &vars.itg_dom[self.dom_id];
        let diff = &vars.scl_dom[self.diff_id];
        let unk = &vars.scl_dom[self.unk_id];

        // iterate over elements
        for eid in 0..dom.num_elem {
            // step 1: assemble local matrix

            // initialize local matrices
            let num_node = dom.elem_node_num[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];
            let gradn_x = &itg.gradn_x[eid];
            let gradn_y = &itg.gradn_y[eid];

            // assemble local matrix
            match num_node {
                3 => {
                    for qid in 0..num_quad {
                        let diff_val = diff.compute_quad(vars, eid, qid, t);
                        let coeff = factor * diff_val * W_TRI3[qid] * jac_det[qid];
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff
                                    * (gradn_x[qid][v] * gradn_x[qid][j]
                                        + gradn_y[qid][v] * gradn_y[qid][j]);
                            }
                        }
                    }
                }
                4 => {
                    for qid in 0..num_quad {
                        let diff_val = diff.compute_quad(vars, eid, qid, t);
                        let coeff = factor * diff_val * W_QUAD4[qid] * jac_det[qid];
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff
                                    * (gradn_x[qid][v] * gradn_x[qid][j]
                                        + gradn_y[qid][v] * gradn_y[qid][j]);
                            }
                        }
                    }
                }
                _ => {
                    panic!("Invalid element type");
                }
            }

            // step 2: assemble global matrix

            // iterate over local matrix entries
            let node_id = &dom.elem_node_id[eid];
            for v in 0..num_node {
                for j in 0..num_node {
                    // get node ids
                    let nid_v = node_id[v];
                    let nid_j = node_id[j];

                    // skip if dirichlet BC
                    if unk.node_dir[nid_v] {
                        continue;
                    }

                    // add to global matrix
                    self.add_a_sclscl(
                        vars,
                        a_triplet,
                        self.unk_id,
                        nid_v,
                        self.drv_id,
                        nid_j,
                        a_loc[v][j],
                    );
                }
            }
        }
    }
}
