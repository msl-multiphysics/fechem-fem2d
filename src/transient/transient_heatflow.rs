use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vec_dom::VectorDomainType;
use crate::base::vec_itf::VectorInterfaceType;
use crate::base::vars::Variables;
use crate::operator::prelude::*;
use crate::transient::transient_base::TransientBase;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct TransientHeatFlow {
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
    pub oper_itr_heat_adv: Vec<(OpSclDomAdvection, OpSclDomSupgSteady, OpSclDomSupgTime)>,
    pub oper_bnd_temp: Vec<OpSclBndDirichlet>,
    pub oper_bnd_hflx: Vec<OpSclBndNeumann>,
    pub oper_bnd_htrn: Vec<OpSclBndTransfer>,
    pub oper_hcnt_itf: Vec<OpSclItfContinuity>,
    pub oper_hres_itf: Vec<OpSclItfTransfer>,

    // momentum transfer fields

    // internal data
    pub itr_fdom: Vec<usize>,            // dom id
    pub itr_vel: HashMap<usize, usize>,  // velocity (unknown)
    pub itr_pres: HashMap<usize, usize>, // pressure (unknown)
    pub itr_den: HashMap<usize, usize>,  // density
    pub itr_visc: HashMap<usize, usize>, // viscosity
    pub itr_fce: HashMap<usize, usize>,  // body force

    // boundary data
    pub vel_bnd: Vec<usize>,              // bnd with velocity BC
    pub vel_vel: HashMap<usize, usize>,   // velocity
    pub pres_bnd: Vec<usize>,             // bnd with pressure BC
    pub pres_pres: HashMap<usize, usize>, // pressure

    // interface data
    pub fcnt_itf: Vec<usize>,                 // itf with continuity BC
    pub fcnt_lmd_vel: HashMap<usize, usize>,  // velocity continuity
    pub fcnt_lmd_pres: HashMap<usize, usize>, // pressure continuity

    // reference pressure
    pub pref_use: bool,
    pub pref_dom: usize,
    pub pref_nid: usize,
    pub pref_pres: f64,

    // operators
    pub oper_itr_flow: Vec<(OpVecDomTime, OpVecDomAdvection, OpVecDomPressure, OpVecDomDiffusion, OpVecDomSource, OpVecDomSupgSteady, OpVecDomSupgTime, OpSclDomDivergence, OpSclDomPspgSteady, OpSclDomPspgTime)>,
    pub oper_bnd_vel: Vec<(OpVecBndDirichlet, OpSclBndDivergence)>,
    pub oper_bnd_pres: Vec<(OpSclBndDirichlet, OpVecBndPressure)>,
    pub oper_fcnt_vel_itf: Vec<OpVecItfContinuity>,
    pub oper_fcnt_pres_itf: Vec<OpSclItfContinuity>,
}

impl TransientHeatFlow {
    pub fn new() -> TransientHeatFlow {
        TransientHeatFlow::default()
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

    pub fn add_flow_dom(&mut self, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) {
        self.itr_fdom.push(dom_id);
        self.itr_vel.insert(dom_id, vel_id);
        self.itr_pres.insert(dom_id, pres_id);
        self.itr_den.insert(dom_id, den_id);
        self.itr_visc.insert(dom_id, visc_id);
        self.itr_fce.insert(dom_id, fce_id);
    }

    pub fn add_vel_bnd(&mut self, bnd_id: usize, vel_id: usize) {
        self.vel_bnd.push(bnd_id);
        self.vel_vel.insert(bnd_id, vel_id);
    }

    pub fn add_pres_bnd(&mut self, bnd_id: usize, pres_id: usize) {
        self.pres_bnd.push(bnd_id);
        self.pres_pres.insert(bnd_id, pres_id);
    }

    pub fn add_fcnt_itf(&mut self, itf_id: usize, lmd_vel_id: usize, lmd_pres_id: usize) {
        self.fcnt_itf.push(itf_id);
        self.fcnt_lmd_vel.insert(itf_id, lmd_vel_id);
        self.fcnt_lmd_pres.insert(itf_id, lmd_pres_id);
    }

    pub fn set_pres_ref(&mut self, dom_id: usize, nid: usize, pres: f64) {
        self.pref_use = true;
        self.pref_dom = dom_id;
        self.pref_nid = nid;
        self.pref_pres = pres;
    }
}

impl TransientBase for TransientHeatFlow {
    fn initial_matrix(&self, vars: &mut Variables) -> usize {
        // assign start indices for unknowns
        let mut xid = 0;

        // flow unknowns
        for (&dom_id, &vel_id) in &self.itr_vel {
            let num_node = vars.dom[dom_id].num_node;
            vars.vec_dom[vel_id].vec_type = VectorDomainType::Unknown { start: xid };
            xid += 2 * num_node;
        }
        for (&dom_id, &pres_id) in &self.itr_pres {
            let num_node = vars.dom[dom_id].num_node;
            vars.scl_dom[pres_id].scl_type = ScalarDomainType::Unknown { start: xid };
            xid += num_node;
        }

        // heat unknowns
        for (&dom_id, &temp_id) in &self.itr_temp {
            let num_node = vars.dom[dom_id].num_node;
            vars.scl_dom[temp_id].scl_type = ScalarDomainType::Unknown { start: xid };
            xid += num_node;
        }

        // flow interface unknowns
        for &itf_id in &self.fcnt_itf {
            let lmd_vel_id = self.fcnt_lmd_vel[&itf_id];
            let num_node = vars.itf[itf_id].num_node;
            vars.vec_itf[lmd_vel_id].vec_type = VectorInterfaceType::Unknown { start: xid };
            xid += 2 * num_node;
        }
        for &itf_id in &self.fcnt_itf {
            let lmd_pres_id = self.fcnt_lmd_pres[&itf_id];
            let num_node = vars.itf[itf_id].num_node;
            vars.scl_itf[lmd_pres_id].scl_type = ScalarInterfaceType::Unknown { start: xid };
            xid += num_node;
        }

        // heat interface unknowns
        for (&itf_id, &lmd_id) in &self.hcnt_lmd {
            let num_node = vars.itf[itf_id].num_node;
            vars.scl_itf[lmd_id].scl_type = ScalarInterfaceType::Unknown { start: xid };
            xid += num_node;
        }

        // return matrix size
        xid
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // list of mesh node ids with dirichlet boundaries
        let mut temp_dir_nid: HashSet<usize> = HashSet::new();
        let mut vel_dir_nid: HashSet<usize> = HashSet::new();
        let mut pres_dir_nid: HashSet<usize> = HashSet::new();

        // iterate over temperature, velocity, and pressure boundaries
        for &bnd_id in &self.temp_bnd {
            for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                temp_dir_nid.insert(mesh_nid);
            }
        }
        for &bnd_id in &self.vel_bnd {
            for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                vel_dir_nid.insert(mesh_nid);
            }
        }
        for &bnd_id in &self.pres_bnd {
            for &mesh_nid in &vars.bnd[bnd_id].node_bnd_mesh_id {
                pres_dir_nid.insert(mesh_nid);
            }
        }
        if self.pref_use {  // reference pressure point
            let dom = &vars.dom[self.pref_dom];
            let mesh_nid = dom.node_dom_mesh_id[self.pref_nid];
            pres_dir_nid.insert(mesh_nid);
        }

        // iterate over heat scalars and flag dirichlet boundaries
        for (&dom_id, &temp_id) in &self.itr_temp {
            let dom = &vars.dom[dom_id];
            let temp = &mut vars.scl_dom[temp_id];
            for dom_nid in 0..dom.num_node {
                let mesh_nid = dom.node_dom_mesh_id[dom_nid];
                temp.node_dir[dom_nid] = temp_dir_nid.contains(&mesh_nid);
            }
        }
        for (&bnd_id, &hflx_id) in &self.hflx_hflx {
            let bnd = &vars.bnd[bnd_id];
            let hflx = &mut vars.scl_bnd[hflx_id];
            for bnd_nid in 0..bnd.num_node {
                let mesh_nid = bnd.node_bnd_mesh_id[bnd_nid];
                hflx.node_dir[bnd_nid] = temp_dir_nid.contains(&mesh_nid);
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

        // iterate over flow fields and flag dirichlet boundaries
        for (&dom_id, &vel_id) in &self.itr_vel {
            let dom = &vars.dom[dom_id];
            let vel = &mut vars.vec_dom[vel_id];
            for dom_nid in 0..dom.num_node {
                let mesh_nid = dom.node_dom_mesh_id[dom_nid];
                vel.node_dir[dom_nid] = vel_dir_nid.contains(&mesh_nid);
            }
        }
        for (&dom_id, &pres_id) in &self.itr_pres {
            let dom = &vars.dom[dom_id];
            let pres = &mut vars.scl_dom[pres_id];
            for dom_nid in 0..dom.num_node {
                let mesh_nid = dom.node_dom_mesh_id[dom_nid];
                pres.node_dir[dom_nid] = pres_dir_nid.contains(&mesh_nid);
            }
        }
        for &itf_id in &self.fcnt_itf {
            let lmd_vel_id = self.fcnt_lmd_vel[&itf_id];
            let itf = &vars.itf[itf_id];
            let lmd = &mut vars.vec_itf[lmd_vel_id];
            for itf_nid in 0..itf.num_node {
                let mesh_nid = itf.node_itf_mesh_id[itf_nid];
                lmd.node_dir[itf_nid] = vel_dir_nid.contains(&mesh_nid);
            }
        }
        for &itf_id in &self.fcnt_itf {
            let lmd_pres_id = self.fcnt_lmd_pres[&itf_id];
            let itf = &vars.itf[itf_id];
            let lmd = &mut vars.scl_itf[lmd_pres_id];
            for itf_nid in 0..itf.num_node {
                let mesh_nid = itf.node_itf_mesh_id[itf_nid];
                lmd.node_dir[itf_nid] = pres_dir_nid.contains(&mesh_nid);
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

            // advective heat transport on overlapping heat+flow domains
            if self.itr_vel.contains_key(&dom_id) {
                let vel_id = self.itr_vel[&dom_id];
                let oper_adv = OpSclDomAdvection::new(dom_id, vlcp_id, vel_id, temp_id);
                let oper_supg = OpSclDomSupgSteady::new(dom_id, vlcp_id, cond_id, vel_id, hsrc_id, temp_id, vec![(cond_id, temp_id)]);
                let oper_supg_time = OpSclDomSupgTime::new(dom_id, vlcp_id, cond_id, vel_id, temp_id);
                self.oper_itr_heat_adv.push((oper_adv, oper_supg, oper_supg_time));
            }
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

        // flow internal operators
        for &dom_id in &self.itr_fdom {
            let vel_id = self.itr_vel[&dom_id];
            let pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let visc_id = self.itr_visc[&dom_id];
            let fce_id = self.itr_fce[&dom_id];

            let oper_time = OpVecDomTime::new(dom_id, den_id, vel_id);
            let oper_adv = OpVecDomAdvection::new(dom_id, den_id, vel_id, vel_id);
            let oper_pres = OpVecDomPressure::new(dom_id, vel_id, pres_id);
            let oper_diff = OpVecDomDiffusion::new(dom_id, visc_id, vel_id, vel_id);
            let oper_src = OpVecDomSource::new(dom_id, fce_id, vel_id);
            let oper_supg = OpVecDomSupgSteady::new(dom_id, den_id, visc_id, vel_id, pres_id, fce_id, vel_id);
            let oper_supg_time = OpVecDomSupgTime::new(dom_id, den_id, visc_id, vel_id, vel_id);
            let oper_div = OpSclDomDivergence::new(dom_id, den_id, vel_id, pres_id);
            let oper_pspg = OpSclDomPspgSteady::new(dom_id, den_id, visc_id, vel_id, pres_id, fce_id, pres_id);
            let oper_pspg_time = OpSclDomPspgTime::new(dom_id, den_id, visc_id, vel_id, pres_id);

            self.oper_itr_flow.push((oper_time, oper_adv, oper_pres, oper_diff, oper_src, oper_supg, oper_supg_time, oper_div, oper_pspg, oper_pspg_time));
        }

        // flow boundary velocity operator
        for &bnd_id in &self.vel_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let bnd_vel_id = self.vel_vel[&bnd_id];
            let oper_dir = OpVecBndDirichlet::new(bnd_id, bnd_vel_id, dom_vel_id);
            let oper_div = OpSclBndDivergence::new(bnd_id, den_id, bnd_vel_id, dom_pres_id);
            self.oper_bnd_vel.push((oper_dir, oper_div));
        }

        // flow boundary pressure operator
        for &bnd_id in &self.pres_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let bnd_pres_id = self.pres_pres[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(bnd_id, bnd_pres_id, dom_pres_id);
            let oper_pres = OpVecBndPressure::new(bnd_id, bnd_pres_id, dom_vel_id);
            self.oper_bnd_pres.push((oper_dir, oper_pres));
        }

        // flow interface continuity operators
        for &itf_id in &self.fcnt_itf {
            let dom1_id = vars.itf[itf_id].dom1_id;
            let dom2_id = vars.itf[itf_id].dom2_id;
            let dom_vel1_id = self.itr_vel[&dom1_id];
            let dom_vel2_id = self.itr_vel[&dom2_id];
            let dom_pres1_id = self.itr_pres[&dom1_id];
            let dom_pres2_id = self.itr_pres[&dom2_id];
            let lmd_vel_id = self.fcnt_lmd_vel[&itf_id];
            let lmd_pres_id = self.fcnt_lmd_pres[&itf_id];
            let oper_cont_vel = OpVecItfContinuity::new(itf_id, lmd_vel_id, dom_vel1_id, dom_vel2_id);
            let oper_cont_pres = OpSclItfContinuity::new(itf_id, lmd_pres_id, dom_pres1_id, dom_pres2_id);
            self.oper_fcnt_vel_itf.push(oper_cont_vel);
            self.oper_fcnt_pres_itf.push(oper_cont_pres);
        }
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, dt: f64) {
        // assemble flow internal data
        for (oper_time, oper_adv, oper_pres, oper_diff, oper_src, oper_supg, oper_supg_time, oper_div, oper_pspg, oper_pspg_time) in &self.oper_itr_flow {
            oper_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
            oper_adv.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_pres.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_diff.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_src.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_supg.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_supg_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
            oper_div.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_pspg.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_pspg_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
        }

        // assemble heat internal data
        for (oper_time, oper_cond, oper_src) in &self.oper_itr_heat {
            oper_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
            oper_cond.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_src.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for (oper_adv, oper_supg, oper_supg_time) in &self.oper_itr_heat_adv {
            oper_adv.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_supg.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_supg_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
        }

        // assemble flow boundary data
        for (oper_dir, oper_div) in &self.oper_bnd_vel {
            oper_dir.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_div.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for (oper_dir, oper_pres) in &self.oper_bnd_pres {
            oper_dir.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_pres.apply(vars, a_triplet, b_vec, t, 1.0);
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

        // assemble flow interface data
        for oper_cont in &self.oper_fcnt_vel_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_cont in &self.oper_fcnt_pres_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble heat interface data
        for oper_cont in &self.oper_hcnt_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_hres in &self.oper_hres_itf {
            oper_hres.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // set reference pressure
        if self.pref_use {
            let pres_id = self.itr_pres[&self.pref_dom];
            let xid = self.pref_nid + match vars.scl_dom[pres_id].scl_type {
                ScalarDomainType::Unknown { start } => start,
                _ => panic!("Expected unknown scalar domain type."),
            };
            a_triplet.push(Triplet::new(xid, xid, 1.0));
            b_vec[xid] = self.pref_pres;
        }
    }
}
