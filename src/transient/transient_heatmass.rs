use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vars::Variables;
use crate::operator::prelude::*;
use crate::transient::transient_base::TransientBase;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct TransientHeatMass {
    // heat transfer fields

    // internal data
    pub itr_hdom: Vec<usize>,            // dom id
    pub itr_temp: HashMap<usize, usize>, // temperature (unknown)
    pub itr_vlcp: HashMap<usize, usize>, // volumetric heat capacity
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
    pub hcnt_itf: Vec<usize>,             // itf with continuity BC
    pub hcnt_lmd: HashMap<usize, usize>,  // continuity
    pub hres_itf: Vec<usize>,             // itf with heat resistance BC
    pub hres_htrn: HashMap<usize, usize>, // heat transfer coefficient

    // operators
    pub oper_itr_heat: Vec<(OpSclDomTime, OpSclDomDiffusion, OpSclDomSource)>,
    pub oper_bnd_temp: Vec<OpSclBndDirichlet>,
    pub oper_bnd_hflx: Vec<OpSclBndNeumann>,
    pub oper_bnd_htrn: Vec<OpSclBndTransfer>,
    pub oper_hcnt_itf: Vec<OpSclItfContinuity>,
    pub oper_hres_itf: Vec<OpSclItfTransfer>,

    // mass transfer fields

    // number of components
    pub num_comp: usize,

    // internal data
    pub itr_mdom: Vec<usize>,             // dom id
    pub itr_conc: Vec<HashMap<usize, usize>>, // concentration (unknown)
    pub itr_diff: Vec<HashMap<usize, HashMap<usize, usize>>>,  // diffusion coefficient
    pub itr_msrc: Vec<HashMap<usize, usize>>, // mass source

    // boundary data
    pub conc_bnd: Vec<Vec<usize>>,             // bnd with concentration BC
    pub conc_conc: Vec<HashMap<usize, usize>>, // concentration
    pub mflx_bnd: Vec<Vec<usize>>,             // bnd with molar flux BC
    pub mflx_mflx: Vec<HashMap<usize, usize>>, // molar flux
    pub mtrn_bnd: Vec<Vec<usize>>,             // bnd with mass transfer BC
    pub mtrn_mtrn: Vec<HashMap<usize, usize>>, // mass transfer coefficient
    pub mtrn_cext: Vec<HashMap<usize, usize>>, // external concentration

    // interface data
    pub mcnt_itf: Vec<Vec<usize>>,             // itf with continuity BC
    pub mcnt_lmd: Vec<HashMap<usize, usize>>,  // continuity
    pub mres_itf: Vec<Vec<usize>>,             // itf with mass resistance BC
    pub mres_mtrn: Vec<HashMap<usize, usize>>, // mass transfer coefficient

    // operators
    pub oper_itr_time: Vec<OpSclDomTimeUnity>,
    pub oper_itr_diff: Vec<OpSclDomDiffusion>,
    pub oper_itr_msrc: Vec<OpSclDomSource>,
    pub oper_bnd_conc: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mflx: Vec<OpSclBndNeumann>,
    pub oper_bnd_mtrn: Vec<OpSclBndTransfer>,
    pub oper_mcnt_itf: Vec<OpSclItfContinuity>,
    pub oper_mres_itf: Vec<OpSclItfTransfer>,
}

impl TransientHeatMass {
    pub fn new(num_comp: usize) -> TransientHeatMass {
        // initialize struct
        let mut phys = TransientHeatMass::default();

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
        phys.mcnt_itf = vec![Vec::new(); num_comp];
        phys.mcnt_lmd = vec![HashMap::new(); num_comp];
        phys.mres_itf = vec![Vec::new(); num_comp];
        phys.mres_mtrn = vec![HashMap::new(); num_comp];

        // return struct
        phys
    }

    pub fn add_heat_dom(&mut self, dom_id: usize, temp_id: usize, vlcp_id: usize, cond_id: usize, hsrc_id: usize) {
        self.itr_hdom.push(dom_id);
        self.itr_temp.insert(dom_id, temp_id);
        self.itr_vlcp.insert(dom_id, vlcp_id);
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

    pub fn add_hcnt_itf(&mut self, itf_id: usize, lmd_id: usize) {
        self.hcnt_itf.push(itf_id);
        self.hcnt_lmd.insert(itf_id, lmd_id);
    }

    pub fn add_hres_itf(&mut self, itf_id: usize, htrn_id: usize) {
        self.hres_itf.push(itf_id);
        self.hres_htrn.insert(itf_id, htrn_id);
    }

    pub fn add_mass_dom(&mut self, comp_i: usize, dom_id: usize, conc_id: usize, diff_ids: HashMap<usize, usize>, msrc_id: usize) {
        if !self.itr_mdom.contains(&dom_id) {
            self.itr_mdom.push(dom_id);
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

    pub fn add_mcnt_itf(&mut self, comp_i: usize, itf_id: usize, lmd_id: usize) {
        self.mcnt_itf[comp_i].push(itf_id);
        self.mcnt_lmd[comp_i].insert(itf_id, lmd_id);
    }

    pub fn add_mres_itf(&mut self, comp_i: usize, itf_id: usize, mtrn_id: usize) {
        self.mres_itf[comp_i].push(itf_id);
        self.mres_mtrn[comp_i].insert(itf_id, mtrn_id);
    }
}

impl TransientBase for TransientHeatMass {
    fn initial_matrix(&self, vars: &mut Variables) -> usize {
        // assign start indices for unknowns
        let mut xid = 0;

        // heat unknowns
        for (&dom_id, &temp_id) in &self.itr_temp {
            let num_node = vars.dom[dom_id].num_node;
            vars.scl_dom[temp_id].scl_type = ScalarDomainType::Unknown { start: xid };
            xid += num_node;
        }

        // mass unknowns
        for comp_i in 0..self.num_comp {
            for (&dom_id, &conc_id) in &self.itr_conc[comp_i] {
                let num_node = vars.dom[dom_id].num_node;
                vars.scl_dom[conc_id].scl_type = ScalarDomainType::Unknown { start: xid };
                xid += num_node;
            }
        }

        // heat interface unknowns
        for (&itf_id, &lmd_id) in &self.hcnt_lmd {
            let num_node = vars.itf[itf_id].num_node;
            vars.scl_itf[lmd_id].scl_type = ScalarInterfaceType::Unknown { start: xid };
            xid += num_node;
        }

        // mass interface unknowns
        for comp_i in 0..self.num_comp {
            for (&itf_id, &lmd_id) in &self.mcnt_lmd[comp_i] {
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
        let mut temp_dir_nid: HashSet<usize> = HashSet::new();

        // iterate over temperature boundaries
        for &bnd_id in &self.temp_bnd {
            for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                temp_dir_nid.insert(mesh_nid);
            }
        }

        // iterate over heat scalars and flag dirichlet boundaries
        // flag per domain so a Dirichlet BC on one domain does not mark the
        // shared interface node as Dirichlet in the neighboring domain
        for (&dom_id, &temp_id) in &self.itr_temp {
            let mut dom_temp_dir: HashSet<usize> = HashSet::new();
            for &bnd_id in &self.temp_bnd {
                if vars.bnd[bnd_id].dom_id == dom_id {
                    for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                        dom_temp_dir.insert(mesh_nid);
                    }
                }
            }
            let dom = &vars.dom[dom_id];
            let temp = &mut vars.scl_dom[temp_id];
            for dom_nid in 0..dom.num_node {
                let mesh_nid = dom.node_dom_mesh_id[dom_nid];
                temp.node_dir[dom_nid] = dom_temp_dir.contains(&mesh_nid);
            }
        }
        for (&bnd_id, &hflx_id) in &self.hflx_hflx {
            let bnd = &vars.bnd[bnd_id];
            let temp_id = self.itr_temp[&bnd.dom_id];
            let temp = &vars.scl_dom[temp_id];
            let hflx = &mut vars.scl_bnd[hflx_id];
            for bnd_nid in 0..bnd.num_node {
                let nid_dom = bnd.node_bnd_dom_id[bnd_nid];
                hflx.node_dir[bnd_nid] = temp.node_dir[nid_dom];
            }
        }
        for (&itf_id, &lmd_id) in &self.hcnt_lmd {
            let itf = &vars.itf[itf_id];
            let lmd = &mut vars.scl_itf[lmd_id];
            for itf_nid in 0..itf.num_node {
                let mesh_nid = itf.node_itf_mesh_id[itf_nid];
                lmd.node_dir[itf_nid] = temp_dir_nid.contains(&mesh_nid);
            }
        }

        // iterate over mass scalars and flag dirichlet boundaries
        for comp_i in 0..self.num_comp {
            let mut conc_dir_nid: HashSet<usize> = HashSet::new();
            for &bnd_id in &self.conc_bnd[comp_i] {
                for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                    conc_dir_nid.insert(mesh_nid);
                }
            }

            for (&dom_id, &conc_id) in &self.itr_conc[comp_i] {
                let mut dom_conc_dir: HashSet<usize> = HashSet::new();
                for &bnd_id in &self.conc_bnd[comp_i] {
                    if vars.bnd[bnd_id].dom_id == dom_id {
                        for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                            dom_conc_dir.insert(mesh_nid);
                        }
                    }
                }
                let dom = &vars.dom[dom_id];
                let conc = &mut vars.scl_dom[conc_id];
                for dom_nid in 0..dom.num_node {
                    let mesh_nid = dom.node_dom_mesh_id[dom_nid];
                    conc.node_dir[dom_nid] = dom_conc_dir.contains(&mesh_nid);
                }
            }
            for (&bnd_id, &mflx_id) in &self.mflx_mflx[comp_i] {
                let bnd = &vars.bnd[bnd_id];
                let conc_id = self.itr_conc[comp_i][&bnd.dom_id];
                let conc = &vars.scl_dom[conc_id];
                let mflx = &mut vars.scl_bnd[mflx_id];
                for bnd_nid in 0..bnd.num_node {
                    let nid_dom = bnd.node_bnd_dom_id[bnd_nid];
                    mflx.node_dir[bnd_nid] = conc.node_dir[nid_dom];
                }
            }
            for (&itf_id, &lmd_id) in &self.mcnt_lmd[comp_i] {
                let itf = &vars.itf[itf_id];
                let lmd = &mut vars.scl_itf[lmd_id];
                for itf_nid in 0..itf.num_node {
                    let mesh_nid = itf.node_itf_mesh_id[itf_nid];
                    lmd.node_dir[itf_nid] = conc_dir_nid.contains(&mesh_nid);
                }
            }
        }
    }

    fn initial_operator(&mut self, vars: &mut Variables) {
        // heat internal operators
        for &dom_id in &self.itr_hdom {
            let temp_id = self.itr_temp[&dom_id];
            let vlcp_id = self.itr_vlcp[&dom_id];
            let cond_id = self.itr_cond[&dom_id];
            let hsrc_id = self.itr_hsrc[&dom_id];
            let oper_time = OpSclDomTime::new(dom_id, vlcp_id, temp_id);
            let oper_cond = OpSclDomDiffusion::new(dom_id, cond_id, temp_id, temp_id);
            let oper_src = OpSclDomSource::new(dom_id, hsrc_id, temp_id);
            self.oper_itr_heat.push((oper_time, oper_cond, oper_src));
        }

        // heat boundary temperature operator
        for &bnd_id in &self.temp_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_temp_id = self.temp_temp[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(bnd_id, bnd_temp_id, dom_temp_id);
            self.oper_bnd_temp.push(oper_dir);
        }

        // heat boundary flux operator
        for &bnd_id in &self.hflx_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_hflx_id = self.hflx_hflx[&bnd_id];
            let oper_neu = OpSclBndNeumann::new(bnd_id, bnd_hflx_id, dom_temp_id);
            self.oper_bnd_hflx.push(oper_neu);
        }

        // heat boundary transfer operator
        for &bnd_id in &self.htrn_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_htrn_id = self.htrn_htrn[&bnd_id];
            let bnd_text_id = self.htrn_text[&bnd_id];
            let oper_htrn = OpSclBndTransfer::new(bnd_id, bnd_htrn_id, bnd_text_id, dom_temp_id);
            self.oper_bnd_htrn.push(oper_htrn);
        }

        // heat interface continuity operator
        for &itf_id in &self.hcnt_itf {
            let dom1_id = vars.itf[itf_id].dom1_id;
            let dom2_id = vars.itf[itf_id].dom2_id;
            let dom_temp1_id = self.itr_temp[&dom1_id];
            let dom_temp2_id = self.itr_temp[&dom2_id];
            let lmd_id = self.hcnt_lmd[&itf_id];
            let oper_cont = OpSclItfContinuity::new(itf_id, lmd_id, dom_temp1_id, dom_temp2_id);
            self.oper_hcnt_itf.push(oper_cont);
        }

        // heat interface heat resistance operator
        for &itf_id in &self.hres_itf {
            let dom1_id = vars.itf[itf_id].dom1_id;
            let dom2_id = vars.itf[itf_id].dom2_id;
            let dom_temp1_id = self.itr_temp[&dom1_id];
            let dom_temp2_id = self.itr_temp[&dom2_id];
            let itf_htrn_id = self.hres_htrn[&itf_id];
            let oper_hres = OpSclItfTransfer::new(itf_id, itf_htrn_id, dom_temp1_id, dom_temp2_id);
            self.oper_hres_itf.push(oper_hres);
        }

        // mass internal operators
        for comp_i in 0..self.num_comp {
            for &dom_id in &self.itr_mdom {
                let conc_id = self.itr_conc[comp_i][&dom_id];
                let oper_time = OpSclDomTimeUnity::new(dom_id, conc_id);
                self.oper_itr_time.push(oper_time);

                let msrc_id = self.itr_msrc[comp_i][&dom_id];
                let oper_src = OpSclDomSource::new(dom_id, msrc_id, conc_id);
                self.oper_itr_msrc.push(oper_src);

                for (&comp_j, &diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let conc_j_id = self.itr_conc[comp_j][&dom_id];
                    let oper_diff = OpSclDomDiffusion::new(dom_id, diff_id, conc_id, conc_j_id);
                    self.oper_itr_diff.push(oper_diff);
                }
            }

            for &bnd_id in &self.conc_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_conc_id = self.conc_conc[comp_i][&bnd_id];
                let oper_dir = OpSclBndDirichlet::new(bnd_id, bnd_conc_id, dom_conc_id);
                self.oper_bnd_conc.push(oper_dir);
            }

            for &bnd_id in &self.mflx_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_mflx_id = self.mflx_mflx[comp_i][&bnd_id];
                let oper_neu = OpSclBndNeumann::new(bnd_id, bnd_mflx_id, dom_conc_id);
                self.oper_bnd_mflx.push(oper_neu);
            }

            for &bnd_id in &self.mtrn_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_mtrn_id = self.mtrn_mtrn[comp_i][&bnd_id];
                let bnd_cext_id = self.mtrn_cext[comp_i][&bnd_id];
                let oper_mtrn = OpSclBndTransfer::new(bnd_id, bnd_mtrn_id, bnd_cext_id, dom_conc_id);
                self.oper_bnd_mtrn.push(oper_mtrn);
            }

            for &itf_id in &self.mcnt_itf[comp_i] {
                let dom1_id = vars.itf[itf_id].dom1_id;
                let dom2_id = vars.itf[itf_id].dom2_id;
                let dom_conc1_id = self.itr_conc[comp_i][&dom1_id];
                let dom_conc2_id = self.itr_conc[comp_i][&dom2_id];
                let lmd_id = self.mcnt_lmd[comp_i][&itf_id];
                let oper_cont = OpSclItfContinuity::new(itf_id, lmd_id, dom_conc1_id, dom_conc2_id);
                self.oper_mcnt_itf.push(oper_cont);
            }

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

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, dt: f64) {
        // assemble heat internal data
        for (oper_time, oper_cond, oper_src) in &self.oper_itr_heat {
            oper_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
            oper_cond.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_src.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble mass internal data
        for oper_time in &self.oper_itr_time {
            oper_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
        }
        for oper_src in &self.oper_itr_msrc {
            oper_src.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_diff in &self.oper_itr_diff {
            oper_diff.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble heat boundary data
        for oper_dir in &self.oper_bnd_temp {
            oper_dir.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_neu in &self.oper_bnd_hflx {
            oper_neu.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_htrn in &self.oper_bnd_htrn {
            oper_htrn.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble mass boundary data
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_neu in &self.oper_bnd_mflx {
            oper_neu.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_mtrn in &self.oper_bnd_mtrn {
            oper_mtrn.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble heat interface data
        for oper_cont in &self.oper_hcnt_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_hres in &self.oper_hres_itf {
            oper_hres.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble mass interface data
        for oper_cont in &self.oper_mcnt_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_mres in &self.oper_mres_itf {
            oper_mres.apply(vars, a_triplet, b_vec, t, 1.0);
        }
    }
}
