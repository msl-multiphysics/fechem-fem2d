use crate::base::vars::Variables;
use crate::operator::oper_base::OperatorBase;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfTransfer {
    // interface
    pub itf_id: usize,

    // scalars
    pub trn_id: usize,  // transfer coefficient
    pub unk1_id: usize, // unknown scalar on domain 1
    pub unk2_id: usize, // unknown scalar on domain 2
}

impl OpSclItfTransfer {
    pub fn new(itf_id: usize, trn_id: usize, unk1_id: usize, unk2_id: usize) -> OpSclItfTransfer {
        // applies a transfer BC at the interface
        // d(m_i * c_i)/dt = -div(N_i) + R_i
        // N1 . n1 = h * (c1 - c2)
        // N2 . n2 = -N1 . n1
        // 
        // trn - transfer coefficient (h)
        // unk1 - unknown scalar on domain 1 (c1)
        // unk2 - unknown scalar on domain 2 (c2)
        // n1 and n2 are the outward normals from each domain
        // N1 is the flux of c1 (due to diffusion, advection, etc.)
        // N2 is the flux of c2 (due to diffusion, advection, etc.)
        
        // create struct
        let mut oper_trn = OpSclItfTransfer::default();
        oper_trn.itf_id = itf_id;
        oper_trn.trn_id = trn_id;
        oper_trn.unk1_id = unk1_id;
        oper_trn.unk2_id = unk2_id;

        // result
        oper_trn
    }
}

impl OperatorBase for OpSclItfTransfer {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>, t: f64, factor: f64) {
        // applies the weak form of the transfer term
        // -(div(N), w)_dom = +(N, grad(w))_dom - (N . n, w)_bnd
        // with N1 . n1 = h * (c1 - c2) and N2 . n2 = -N1 . n1
        //
        // let A (in Ax = b) be the RHS of the PDE and b in the LHS
        // c1: add -(h * (c1 - c2), w)_itf to A
        // c2: add +(h * (c1 - c2), w)_itf to A

        // get objects
        let itf = &vars.itf[self.itf_id];
        let itg = &vars.itg_itf[self.itf_id];
        let trn_scl = &vars.scl_itf[self.trn_id];
        let unk1 = &vars.scl_dom[self.unk1_id];
        let unk2 = &vars.scl_dom[self.unk2_id];

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
                let trn = trn_scl.compute_quad(eid, qid, t);
                let coeff = -factor * quad_w[qid] * trn * jac_det[qid].sqrt();
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
                let nid_itf_v = node_id[v];
                let nid1_v = itf.node_itf_dom1_id[nid_itf_v];
                let nid2_v = itf.node_itf_dom2_id[nid_itf_v];

                // add to global matrix
                for j in 0..num_node {
                    // get column indices
                    let nid_itf_j = node_id[j];
                    let nid1_j = itf.node_itf_dom1_id[nid_itf_j];
                    let nid2_j = itf.node_itf_dom2_id[nid_itf_j];

                    // domain 1 rows
                    if !unk1.node_dir[nid1_v] {
                        self.add_a_scldom(vars, a_triplet, self.unk1_id, nid1_v, self.unk1_id, nid1_j, a_loc[v][j]);
                        self.add_a_scldom(vars, a_triplet, self.unk1_id, nid1_v, self.unk2_id, nid2_j, -a_loc[v][j]);
                    }

                    // domain 2 rows
                    if !unk2.node_dir[nid2_v] {
                        self.add_a_scldom(vars, a_triplet, self.unk2_id, nid2_v, self.unk1_id, nid1_j, -a_loc[v][j]);
                        self.add_a_scldom(vars, a_triplet, self.unk2_id, nid2_v, self.unk2_id, nid2_j, a_loc[v][j]);
                    }
                }
            }

        }
    }
}
