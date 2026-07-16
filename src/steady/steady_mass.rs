use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vars::Variables;
use crate::operator::prelude::*;
use crate::steady::steady_base::SteadyBase;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct SteadyMass {
    // number of components
    pub num_comp: usize,

    // internal data
    // concentration (c_i): itr_conc[comp_i][dom_id] -> conc_id
    // diffusion coefficient (D_ij): itr_diff[comp_i][dom_id][comp_j] -> diff_id
    // reaction rate (R_i): itr_msrc[comp_i][dom_id] -> msrc_id
    pub itr_dom: Vec<usize>,             // dom id
    pub itr_conc: Vec<HashMap<usize, usize>>, // concentration (unknown)
    pub itr_diff: Vec<HashMap<usize, HashMap<usize, usize>>>,  // diffusion coefficient
    pub itr_msrc: Vec<HashMap<usize, usize>>, // mass source

    // boundary data
    // concentration (c_i): conc_conc[comp_i][dom1d_id] -> scl1d_conc_id
    // molar flux (N_i): mflx_mflx[comp_i][dom1d_id] -> scl1d_mflx_id
    pub conc_bnd: Vec<Vec<usize>>,             // bnd with concentration BC
    pub conc_conc: Vec<HashMap<usize, usize>>, // concentration
    pub mflx_bnd: Vec<Vec<usize>>,             // bnd with molar flux BC
    pub mflx_mflx: Vec<HashMap<usize, usize>>, // molar flux
    pub mtrn_bnd: Vec<Vec<usize>>,             // bnd with mass transfer BC
    pub mtrn_mtrn: Vec<HashMap<usize, usize>>, // mass transfer coefficient
    pub mtrn_cext: Vec<HashMap<usize, usize>>, // external concentration

    // interface data
    pub cont_itf: Vec<Vec<usize>>,             // itf with continuity BC
    pub cont_lmd: Vec<HashMap<usize, usize>>,  // continuity
    pub mres_itf: Vec<Vec<usize>>,             // itf with mass resistance BC
    pub mres_mtrn: Vec<HashMap<usize, usize>>, // mass transfer coefficient

    // operators
    pub oper_itr_diff: Vec<OpSclDomDiffusion>,
    pub oper_itr_msrc: Vec<OpSclDomSource>,
    pub oper_bnd_conc: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mflx: Vec<OpSclBndNeumann>,
    pub oper_bnd_mtrn: Vec<OpSclBndTransfer>,
    pub oper_cont_itf: Vec<OpSclItfContinuity>,
    pub oper_mres_itf: Vec<OpSclItfTransfer>,
}

impl SteadyMass {
    pub fn new(num_comp: usize) -> SteadyMass {
        // initialize struct
        let mut phys = SteadyMass::default();

        // resize vectors
        phys.num_comp = num_comp;
        phys.itr_conc = vec![HashMap::new(); num_comp];
        phys.itr_diff = vec![HashMap::new(); num_comp];
        phys.itr_msrc = vec![HashMap::new(); num_comp];
        phys.conc_bnd = vec![Vec::new(); num_comp];
        phys.conc_conc = vec![HashMap::new(); num_comp];
        phys.mflx_bnd = vec![Vec::new(); num_comp];
        phys.mflx_mflx = vec![HashMap::new(); num_comp];
        phys.mtrn_bnd = vec![Vec::new(); num_comp];
        phys.mtrn_mtrn = vec![HashMap::new(); num_comp];
        phys.mtrn_cext = vec![HashMap::new(); num_comp];
        phys.cont_itf = vec![Vec::new(); num_comp];
        phys.cont_lmd = vec![HashMap::new(); num_comp];
        phys.mres_itf = vec![Vec::new(); num_comp];
        phys.mres_mtrn = vec![HashMap::new(); num_comp];

        // return struct
        phys
    }

    pub fn add_mass_dom(&mut self, comp_i: usize, dom_id: usize, conc_id: usize, diff_ids: HashMap<usize, usize>, msrc_id: usize) {
        if !self.itr_dom.contains(&dom_id) {
            self.itr_dom.push(dom_id);
        }
        self.itr_conc[comp_i].insert(dom_id, conc_id);
        self.itr_diff[comp_i].insert(dom_id, diff_ids);
        self.itr_msrc[comp_i].insert(dom_id, msrc_id);
    }

    pub fn add_conc_bnd(&mut self, comp_i: usize, bnd_id: usize, conc_id: usize) {
        self.conc_bnd[comp_i].push(bnd_id);
        self.conc_conc[comp_i].insert(bnd_id, conc_id);
    }

    pub fn add_mflx_bnd(&mut self, comp_i: usize, bnd_id: usize, mflx_id: usize) {
        self.mflx_bnd[comp_i].push(bnd_id);
        self.mflx_mflx[comp_i].insert(bnd_id, mflx_id);
    }

    pub fn add_mtrn_bnd(&mut self, comp_i: usize, bnd_id: usize, mtrn_id: usize, cext_id: usize) {
        self.mtrn_bnd[comp_i].push(bnd_id);
        self.mtrn_mtrn[comp_i].insert(bnd_id, mtrn_id);
        self.mtrn_cext[comp_i].insert(bnd_id, cext_id);
    }

    pub fn add_cont_itf(&mut self, comp_i: usize, itf_id: usize, lmd_id: usize) {
        self.cont_itf[comp_i].push(itf_id);
        self.cont_lmd[comp_i].insert(itf_id, lmd_id);
    }

    pub fn add_mres_itf(&mut self, comp_i: usize, itf_id: usize, mtrn_id: usize) {
        self.mres_itf[comp_i].push(itf_id);
        self.mres_mtrn[comp_i].insert(itf_id, mtrn_id);
    }
}

impl SteadyBase for SteadyMass {
    fn initial_matrix(&self, vars: &mut Variables) -> usize {
        // assign start indices for unknowns
        let mut xid = 0;
        for comp_i in 0..self.num_comp {
            for (&dom_id, &conc_id) in &self.itr_conc[comp_i] {
                let num_node = vars.dom[dom_id].num_node;
                vars.scl_dom[conc_id].scl_type = ScalarDomainType::Unknown { start: xid };
                xid += num_node;
            }
            for (&itf_id, &lmd_id) in &self.cont_lmd[comp_i] {
                let num_node = vars.itf[itf_id].num_node;
                vars.scl_itf[lmd_id].scl_type = ScalarInterfaceType::Unknown { start: xid };
                xid += num_node;
            }
        }

        // return matrix size
        xid
    }
    
    fn initial_dirichlet(&self, vars: &mut Variables) {
        // list of mesh node ids with dirichlet boundaries
        let mut dir_nid: HashSet<usize> = HashSet::new();

        // iterate over concentration boundaries
        for comp_i in 0..self.num_comp {
            for &bnd_id in &self.conc_bnd[comp_i] {
                for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                    dir_nid.insert(mesh_nid);
                }
            }
        }

        // iterate over scalars and flag dirichlet boundaries
        for comp_i in 0..self.num_comp {
            for (&dom_id, &conc_id) in &self.itr_conc[comp_i] {
                let dom = &vars.dom[dom_id];
                let conc = &mut vars.scl_dom[conc_id];
                for dom_nid in 0..dom.num_node {
                    let mesh_nid = dom.node_dom_mesh_id[dom_nid];
                    conc.node_dir[dom_nid] = dir_nid.contains(&mesh_nid);
                }
            }
            for (&bnd_id, &mflx_id) in &self.mflx_mflx[comp_i] {
                let bnd = &vars.bnd[bnd_id];
                let mflx = &mut vars.scl_bnd[mflx_id];
                for bnd_nid in 0..bnd.num_node {
                    let mesh_nid = bnd.node_bnd_mesh_id[bnd_nid];
                    mflx.node_dir[bnd_nid] = dir_nid.contains(&mesh_nid);
                }
            }
            for (&itf_id, &lmd_id) in &self.cont_lmd[comp_i] {
                let itf = &vars.itf[itf_id];
                let lmd = &mut vars.scl_itf[lmd_id];
                for itf_nid in 0..itf.num_node {
                    let mesh_nid = itf.node_itf_mesh_id[itf_nid];
                    lmd.node_dir[itf_nid] = dir_nid.contains(&mesh_nid);
                }
            }
        }

    }
    
    fn initial_operator(&mut self, vars: &mut Variables) {
        // iterate over components
        for comp_i in 0..self.num_comp {
            // internal operator
            for &dom_id in &self.itr_dom {
                // source operator -> comp_i
                let conc_id = self.itr_conc[comp_i][&dom_id];
                let msrc_id = self.itr_msrc[comp_i][&dom_id];
                let oper_src = OpSclDomSource::new(dom_id, msrc_id, conc_id);
                self.oper_itr_msrc.push(oper_src);

                // diffusion operator -> (comp_i, comp_j)
                for (&comp_j, &diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let conc_j_id = self.itr_conc[comp_j][&dom_id];
                    let oper_diff = OpSclDomDiffusion::new(dom_id, diff_id, conc_id, conc_j_id);
                    self.oper_itr_diff.push(oper_diff);
                }
            }

            // boundary concentration operator
            for &bnd_id in &self.conc_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_conc_id = self.conc_conc[comp_i][&bnd_id];
                let oper_dir = OpSclBndDirichlet::new(bnd_id, bnd_conc_id, dom_conc_id);
                self.oper_bnd_conc.push(oper_dir);
            }

            // boundary flux operator
            for &bnd_id in &self.mflx_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_mflx_id = self.mflx_mflx[comp_i][&bnd_id];
                let oper_neu = OpSclBndNeumann::new(bnd_id, bnd_mflx_id, dom_conc_id);
                self.oper_bnd_mflx.push(oper_neu);
            }

            // boundary transfer operator
            for &bnd_id in &self.mtrn_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_mtrn_id = self.mtrn_mtrn[comp_i][&bnd_id];
                let bnd_cext_id = self.mtrn_cext[comp_i][&bnd_id];
                let oper_mtrn = OpSclBndTransfer::new(bnd_id, bnd_mtrn_id, bnd_cext_id, dom_conc_id);
                self.oper_bnd_mtrn.push(oper_mtrn);
            }

            // interface continuity operator
            for &itf_id in &self.cont_itf[comp_i] {
                let dom1_id = vars.itf[itf_id].dom1_id;
                let dom2_id = vars.itf[itf_id].dom2_id;
                let dom_conc1_id = self.itr_conc[comp_i][&dom1_id];
                let dom_conc2_id = self.itr_conc[comp_i][&dom2_id];
                let lmd_id = self.cont_lmd[comp_i][&itf_id];
                let oper_cont = OpSclItfContinuity::new(itf_id, lmd_id, dom_conc1_id, dom_conc2_id);
                self.oper_cont_itf.push(oper_cont);
            }

            // interface mass resistance operator
            for &itf_id in &self.mres_itf[comp_i] {
                let dom1_id = vars.itf[itf_id].dom1_id;
                let dom2_id = vars.itf[itf_id].dom2_id;
                let dom_conc1_id = self.itr_conc[comp_i][&dom1_id];
                let dom_conc2_id = self.itr_conc[comp_i][&dom2_id];
                let itf_mtrn_id = self.mres_mtrn[comp_i][&itf_id];
                let oper_mres = OpSclItfTransfer::new(itf_id, itf_mtrn_id, dom_conc1_id, dom_conc2_id);
                self.oper_mres_itf.push(oper_mres);
            }
        }
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) {
        // assemble internal data
        for oper_src in &self.oper_itr_msrc {
            oper_src.apply(vars, a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_diff in &self.oper_itr_diff {
            oper_diff.apply(vars, a_triplet, b_vec, 0.0, 1.0);
        }

        // assemble boundary data
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply(vars, a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_neu in &self.oper_bnd_mflx {
            oper_neu.apply(vars, a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_mtrn in &self.oper_bnd_mtrn {
            oper_mtrn.apply(vars, a_triplet, b_vec, 0.0, 1.0);
        }

        // assemble interface data
        for oper_cont in &self.oper_cont_itf {
            oper_cont.apply(vars, a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_mres in &self.oper_mres_itf {
            oper_mres.apply(vars, a_triplet, b_vec, 0.0, 1.0);
        }

    }
}
