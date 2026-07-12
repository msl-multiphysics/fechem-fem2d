use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecItfContinuity {
    // interface
    pub itf_id: usize,

    // vectors
    pub lmd_id: usize,  // lagrange multiplier
    pub unk1_id: usize, // unknown vector on domain 1
    pub unk2_id: usize, // unknown vector on domain 2
}

impl OpVecItfContinuity {
    pub fn new(itf_id: usize, lmd_id: usize, unk1_id: usize, unk2_id: usize) -> OpVecItfContinuity {
        // applies value and flux continuity at the interface
        // v1 = v2
        // T1 . n = T2 . n
        // 
        // lmd - lagrange multiplier
        // unk1 - unknown vector on domain 1 (v1)
        // unk2 - unknown vector on domain 2 (v2)
        // n is the outward normal from each interface
        // T1 is the flux of v1 (due to diffusion, advection, etc.)
        // T2 is the flux of v2 (due to diffusion, advection, etc.)
        
        // create struct
        let mut oper_cont = OpVecItfContinuity::default();
        oper_cont.itf_id = itf_id;
        oper_cont.lmd_id = lmd_id;
        oper_cont.unk1_id = unk1_id;
        oper_cont.unk2_id = unk2_id;

        // result
        oper_cont
    }
}

impl OperatorBase for OpVecItfContinuity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, _t: f64, factor: f64) {
        // applies the continuity constraint
        // ((v1 - v2), l)_itf
        // v1 and v2 are approximated with basis functions w_v
        // l is approximated with basis functions w_l
        // 
        // compute variations wrt each domain
        // v1: +(w_v, l)_itf
        // v2: -(w_v, l)_itf
        // l: +((v1 - v2), w_l)_itf
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // v1: add +(w_v, l)_itf to A
        // v2: add -(w_v, l)_itf to A
        // l: add +((v1 - v2), w_l)_itf to A

        // get objects
        let itf = &vars.itf[self.itf_id];
        let itg = &vars.itg_itf[self.itf_id];
        let unk1 = &vars.vec_dom[self.unk1_id];
        let unk2 = &vars.vec_dom[self.unk2_id];

        // iterate over elements
        for eid in 0..itf.num_elem {
            // step 1: assemble local matrix

            // initialize local matrix
            let num_node = itf.elem_node[eid];
            let mut a_loc = vec![vec![0.0; num_node]; num_node];

            // get quadrature point data
            let num_quad = itg.num_quad[eid];
            let quad_w = &itg.quad_w[eid];
            let quad_n = &itg.quad_n[eid];
            let jac_det = &itg.jac_det[eid];

            // assemble local matrix
            for qid in 0..num_quad {
                let coeff = factor * quad_w[qid] * jac_det[qid].sqrt();
                for v in 0..num_node {
                    for j in 0..num_node {
                        a_loc[v][j] += coeff * quad_n[qid][v] * quad_n[qid][j];
                    }
                }
            }

            // step 2: add to global matrix

            // iterate over nodes in element
            let node_id = &itf.elem_node1_id[eid];
            for v in 0..num_node {
                // get row indices
                let nid_lmd_v = node_id[v];
                let nid1_v = itf.node_itf_dom1_id[nid_lmd_v];
                let nid2_v = itf.node_itf_dom2_id[nid_lmd_v];

                // flag dirichlet boundaries
                let add_unk1_lmd = !unk1.node_dir[nid1_v] || unk2.node_dir[nid2_v];
                let add_unk2_lmd = !unk2.node_dir[nid2_v] || unk1.node_dir[nid1_v];

                // add to global matrix
                for j in 0..num_node {
                    // get column indices
                    let nid_lmd_j = node_id[j];
                    let nid1_j = itf.node_itf_dom1_id[nid_lmd_j];
                    let nid2_j = itf.node_itf_dom2_id[nid_lmd_j];

                    // continuity constraint (x and y components)
                    for comp in 0..2 {
                        self.add_a_vecitf_vecdom(vars, a_triplet, self.lmd_id, comp, nid_lmd_v, self.unk1_id, comp, nid1_j, a_loc[v][j]);
                        self.add_a_vecitf_vecdom(vars, a_triplet, self.lmd_id, comp, nid_lmd_v, self.unk2_id, comp, nid2_j, -a_loc[v][j]);

                        // transpose terms for lagrange multipliers
                        if add_unk1_lmd {
                            self.add_a_vecdom_vecitf(vars, a_triplet, self.unk1_id, comp, nid1_v, self.lmd_id, comp, nid_lmd_j, a_loc[v][j]);
                        }
                        if add_unk2_lmd {
                            self.add_a_vecdom_vecitf(vars, a_triplet, self.unk2_id, comp, nid2_v, self.lmd_id, comp, nid_lmd_j, -a_loc[v][j]);
                        }
                    }
                }
            }

        }
    }
}
