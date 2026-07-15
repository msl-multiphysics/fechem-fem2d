use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vars::Variables;
use crate::operator::prelude::*;
use crate::steady::steady_base::SteadyBase;
use faer::Col;
use faer::sparse::{SparseColMat, Triplet};
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct SteadyHeat {
    // internal data
    pub itr_dom: Vec<usize>,             // dom id
    pub itr_temp: HashMap<usize, usize>, // temperature (unknown)
    pub itr_cond: HashMap<usize, usize>, // thermal conductivity
    pub itr_hsrc: HashMap<usize, usize>, // heat source

    // boundary data
    pub temp_bnd: Vec<usize>,             // bnd with temperature BC
    pub temp_temp: HashMap<usize, usize>, // temperature
    pub hflx_bnd: Vec<usize>,             // bnd with heat flux BC
    pub hflx_hflx: HashMap<usize, usize>, // heat flux
    pub htrn_bnd: Vec<usize>,             // bnd with heat transfer BC
    pub htrn_htrn: HashMap<usize, usize>, // heat transfer coefficient
    pub htrn_text: HashMap<usize, usize>, // external temperature

    // interface data
    pub cont_itf: Vec<usize>,             // itf with continuity BC
    pub cont_lmd: HashMap<usize, usize>,  // continuity
    pub hres_itf: Vec<usize>,             // itf with heat resistance BC
    pub hres_htrn: HashMap<usize, usize>, // heat transfer coefficient

    // operators
    pub oper_itr: Vec<(OpSclDomDiffusion, OpSclDomSource)>,
    pub oper_bnd_temp: Vec<OpSclBndDirichlet>,
    pub oper_bnd_hflx: Vec<OpSclBndNeumann>,
    pub oper_bnd_htrn: Vec<OpSclBndTransfer>,
    pub oper_cont_itf: Vec<OpSclItfContinuity>,
    pub oper_hres_itf: Vec<OpSclItfTransfer>,
}

impl SteadyHeat {
    pub fn new() -> SteadyHeat {
        SteadyHeat::default()
    }

    pub fn add_heat_dom(&mut self, dom_id: usize, temp_id: usize, cond_id: usize, hsrc_id: usize) {
        self.itr_dom.push(dom_id);
        self.itr_temp.insert(dom_id, temp_id);
        self.itr_cond.insert(dom_id, cond_id);
        self.itr_hsrc.insert(dom_id, hsrc_id);
    }

    pub fn add_temp_bnd(&mut self, bnd_id: usize, temp_id: usize) {
        self.temp_bnd.push(bnd_id);
        self.temp_temp.insert(bnd_id, temp_id);
    }

    pub fn add_hflx_bnd(&mut self, bnd_id: usize, hflx_id: usize) {
        self.hflx_bnd.push(bnd_id);
        self.hflx_hflx.insert(bnd_id, hflx_id);
    }

    pub fn add_htrn_bnd(&mut self, bnd_id: usize, htrn_id: usize, text_id: usize) {
        self.htrn_bnd.push(bnd_id);
        self.htrn_htrn.insert(bnd_id, htrn_id);
        self.htrn_text.insert(bnd_id, text_id);
    }

    pub fn add_cont_itf(&mut self, itf_id: usize, lmd_id: usize) {
        self.cont_itf.push(itf_id);
        self.cont_lmd.insert(itf_id, lmd_id);
    }

    pub fn add_hres_itf(&mut self, itf_id: usize, htrn_id: usize) {
        self.hres_itf.push(itf_id);
        self.hres_htrn.insert(itf_id, htrn_id);
    }
}

impl SteadyBase for SteadyHeat {
    fn initial_matrix(&self, vars: &mut Variables) -> usize {
        // assign start indices for unknowns
        let mut xid = 0;
        for (&dom_id, &temp_id) in &self.itr_temp {
            let num_node = vars.dom[dom_id].num_node;
            vars.scl_dom[temp_id].scl_type = ScalarDomainType::Unknown { start: xid };
            xid += num_node;
        }
        for (&itf_id, &lmd_id) in &self.cont_lmd {
            let num_node = vars.itf[itf_id].num_node;
            vars.scl_itf[lmd_id].scl_type = ScalarInterfaceType::Unknown { start: xid };
            xid += num_node;
        }

        // return matrix size
        xid
    }
    
    fn initial_dirichlet(&self, vars: &mut Variables) {
        // list of mesh node ids with dirichlet boundaries
        let mut dir_nid: HashSet<usize> = HashSet::new();

        // iterate over temperature boundaries
        for &bnd_id in &self.temp_bnd {
            for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                dir_nid.insert(mesh_nid);
            }
        }

        // iterate over scalars and flag dirichlet boundaries
        for (&dom_id, &temp_id) in &self.itr_temp {
            let dom = &vars.dom[dom_id];
            let temp = &mut vars.scl_dom[temp_id];
            for dom_nid in 0..dom.num_node {
                let mesh_nid = dom.node_dom_mesh_id[dom_nid];
                temp.node_dir[dom_nid] = dir_nid.contains(&mesh_nid);
            }
        }
        for (&bnd_id, &hflx_id) in &self.hflx_hflx {
            let bnd = &vars.bnd[bnd_id];
            let hflx = &mut vars.scl_bnd[hflx_id];
            for bnd_nid in 0..bnd.num_node {
                let mesh_nid = bnd.node_bnd_mesh_id[bnd_nid];
                hflx.node_dir[bnd_nid] = dir_nid.contains(&mesh_nid);
            }
        }
        for (&itf_id, &lmd_id) in &self.cont_lmd {
            let itf = &vars.itf[itf_id];
            let lmd = &mut vars.scl_itf[lmd_id];
            for itf_nid in 0..itf.num_node {
                let mesh_nid = itf.node_itf_mesh_id[itf_nid];
                lmd.node_dir[itf_nid] = dir_nid.contains(&mesh_nid);
            }
        }
    }
    
    fn initial_operator(&mut self, vars: &mut Variables) {
        // internal operator
        for &dom_id in &self.itr_dom {
            let temp_id = self.itr_temp[&dom_id];
            let cond_id = self.itr_cond[&dom_id];
            let hsrc_id = self.itr_hsrc[&dom_id];
            let oper_cond = OpSclDomDiffusion::new(dom_id, cond_id, temp_id, temp_id);
            let oper_src = OpSclDomSource::new(dom_id, hsrc_id, temp_id);
            self.oper_itr.push((oper_cond, oper_src));
        }

        // boundary temperature operator
        for &bnd_id in &self.temp_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_temp_id = self.temp_temp[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(bnd_id, bnd_temp_id, dom_temp_id);
            self.oper_bnd_temp.push(oper_dir);
        }

        // boundary flux operator
        for &bnd_id in &self.hflx_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_hflx_id = self.hflx_hflx[&bnd_id];
            let oper_neu = OpSclBndNeumann::new(bnd_id, bnd_hflx_id, dom_temp_id);
            self.oper_bnd_hflx.push(oper_neu);
        }

        // boundary transfer operator
        for &bnd_id in &self.htrn_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_htrn_id = self.htrn_htrn[&bnd_id];
            let bnd_text_id = self.htrn_text[&bnd_id];
            let oper_htrn = OpSclBndTransfer::new(bnd_id, bnd_htrn_id, bnd_text_id, dom_temp_id);
            self.oper_bnd_htrn.push(oper_htrn);
        }

        // interface continuity operator
        for &itf_id in &self.cont_itf {
            let dom1_id = vars.itf[itf_id].dom1_id;
            let dom2_id = vars.itf[itf_id].dom2_id;
            let dom_temp1_id = self.itr_temp[&dom1_id];
            let dom_temp2_id = self.itr_temp[&dom2_id];
            let lmd_id = self.cont_lmd[&itf_id];
            let oper_cont = OpSclItfContinuity::new(itf_id, lmd_id, dom_temp1_id, dom_temp2_id);
            self.oper_cont_itf.push(oper_cont);
        }

        // interface heat resistance operator
        for &itf_id in &self.hres_itf {
            let dom1_id = vars.itf[itf_id].dom1_id;
            let dom2_id = vars.itf[itf_id].dom2_id;
            let dom_temp1_id = self.itr_temp[&dom1_id];
            let dom_temp2_id = self.itr_temp[&dom2_id];
            let itf_htrn_id = self.hres_htrn[&itf_id];
            let oper_hres = OpSclItfTransfer::new(itf_id, itf_htrn_id, dom_temp1_id, dom_temp2_id);
            self.oper_hres_itf.push(oper_hres);
        }
    }

    fn assemble_matrix(&self, vars: &Variables, a_mat: &mut SparseColMat<usize, f64>, b_vec: &mut Col<f64>, mat_size: usize) {
        // initialize triplet for matrix assembly
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        *b_vec = Col::zeros(mat_size);

        // assemble internal data
        for (oper_cond, oper_src) in &self.oper_itr {
            oper_cond.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
            oper_src.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }

        // assemble boundary data
        for oper_dir in &self.oper_bnd_temp {
            oper_dir.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_neu in &self.oper_bnd_hflx {
            oper_neu.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_htrn in &self.oper_bnd_htrn {
            oper_htrn.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }

        // assemble interface data
        for oper_cont in &self.oper_cont_itf {
            oper_cont.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_hres in &self.oper_hres_itf {
            oper_hres.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }

        // create sparse matrix from triplet
        *a_mat = SparseColMat::try_new_from_triplets(mat_size, mat_size, &a_triplet).expect("Failed to create sparse matrix from triplets.");
    }
}
