use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use crate::shape::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfContinuity {
    // interface
    pub itf_id: usize,

    // scalars
    pub lmd_id: usize,  // lagrange multiplier
    pub unk1_id: usize, // unknown scalar on domain 1
    pub unk2_id: usize, // unknown scalar on domain 2
}

impl OpSclItfContinuity {
    pub fn new(itf_id: usize, lmd_id: usize, unk1_id: usize, unk2_id: usize) -> OpSclItfContinuity {
        // create struct
        let mut oper_cont = OpSclItfContinuity::default();
        oper_cont.itf_id = itf_id;
        oper_cont.lmd_id = lmd_id;
        oper_cont.unk1_id = unk1_id;
        oper_cont.unk2_id = unk2_id;

        // result
        oper_cont
    }
}

impl OperatorBase for OpSclItfContinuity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // get objects
        let itf = &vars.itf[self.itf_id];
        let itg = &vars.itg_itf[self.itf_id];
        let unk1 = &vars.scl_dom[self.unk1_id];
        let unk2 = &vars.scl_dom[self.unk2_id];

        // iterate over elements
        for eid in 0..itf.num_elem {
            // step 1: assemble local matrix

            // initialize local matrix
            let num_node = itf.elem_node_num[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];

            // get integral data
            let num_quad = itg.num_quad[eid];
            let jac_det = &itg.jac_det[eid];

            // assemble local matrix
            match num_node {
                2 => {
                    for qid in 0..num_quad {
                        let n = lin2_eval(A_LIN2[qid]);
                        let coeff = factor * W_LIN2[qid] * jac_det[qid].sqrt();
                        for v in 0..num_node {
                            for j in 0..num_node {
                                a_loc[v][j] += coeff * n[v] * n[j];
                            }
                        }
                    }
                }
                _ => {
                    panic!("Invalid element type");
                }
            }

            // step 2: assemble global matrix

            // iterate over rows
            let node_id = &itf.elem_node1_id[eid];
            for v in 0..num_node {
                // get row indices
                let nid_lmd_v = node_id[v];
                let nid1_v = itf.node_itf_dom1_id[nid_lmd_v];
                let nid2_v = itf.node_itf_dom2_id[nid_lmd_v];

                // flag dirichlet boundaries
                let add_unk1_lmd = !unk1.node_dir[nid1_v] || unk2.node_dir[nid2_v];
                let add_unk2_lmd = !unk2.node_dir[nid2_v] || unk1.node_dir[nid1_v];

                // iterate over columns
                for j in 0..num_node {
                    // get column indices
                    let nid_lmd_j = node_id[j];
                    let nid1_j = itf.node_itf_dom1_id[nid_lmd_j];
                    let nid2_j = itf.node_itf_dom2_id[nid_lmd_j];

                    // get matrix entry
                    let a_vj = a_loc[v][j];

                    // continuity constraint on PDEs
                    self.add_a_sclitf_scldom(vars, a_triplet, self.lmd_id, nid_lmd_v, self.unk1_id, nid1_j, a_vj);
                    self.add_a_sclitf_scldom(vars, a_triplet, self.lmd_id, nid_lmd_v, self.unk2_id, nid2_j, -a_vj);

                    // transpose terms for lagrange multipliers
                    if add_unk1_lmd {
                        self.add_a_scldom_sclitf(vars, a_triplet, self.unk1_id, nid1_v, self.lmd_id, nid_lmd_j, a_vj);
                    }
                    if add_unk2_lmd {
                        self.add_a_scldom_sclitf(vars, a_triplet, self.unk2_id, nid2_v, self.lmd_id, nid_lmd_j, -a_vj);
                    }
                }
            }
        }
    }
}
