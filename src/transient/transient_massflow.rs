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
pub struct TransientMassFlow {
    // mass transfer fields

    // number of components
    pub num_comp: usize,

    // internal data
    // concentration (c_i): itr_conc[comp_i][dom_id] -> conc_id
    // diffusion coefficient (D_ij): itr_diff[comp_i][dom_id][comp_j] -> diff_id
    // reaction rate (R_i): itr_msrc[comp_i][dom_id] -> msrc_id
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
    pub mout_bnd: Vec<Vec<usize>>,             // bnd with mass outflow BC

    // interface data
    pub mcnt_itf: Vec<Vec<usize>>,             // itf with continuity BC
    pub mcnt_lmd: Vec<HashMap<usize, usize>>,  // continuity
    pub mres_itf: Vec<Vec<usize>>,             // itf with mass resistance BC
    pub mres_mtrn: Vec<HashMap<usize, usize>>, // mass transfer coefficient

    // operators
    pub oper_itr_time: Vec<OpSclDomTimeUnity>,
    pub oper_itr_diff: Vec<OpSclDomDiffusion>,
    pub oper_itr_msrc: Vec<OpSclDomSource>,
    pub oper_itr_mass_adv: Vec<(OpSclDomAdvectionUnity, OpSclDomSupgSteadyUnity, OpSclDomSupgTimeUnity)>,
    pub oper_bnd_conc: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mflx: Vec<OpSclBndNeumann>,
    pub oper_bnd_mtrn: Vec<OpSclBndTransfer>,
    pub oper_bnd_mout: Vec<OpSclBndOutflowUnity>,
    pub oper_mcnt_itf: Vec<OpSclItfContinuity>,
    pub oper_mres_itf: Vec<OpSclItfTransfer>,

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
    pub oper_itr_flow: Vec<(OpVecDomTime, OpVecDomAdvection, OpVecDomPressure, OpVecDomDiffusion, OpVecDomSource, OpVecDomSupgSteady, OpVecDomSupgTime, OpSclDomDensityTime, OpSclDomDivergence, OpSclDomPspgSteady, OpSclDomPspgTime)>,
    pub oper_bnd_vel: Vec<(OpVecBndDirichlet, OpSclBndDivergence)>,
    pub oper_bnd_pres: Vec<(OpSclBndDirichlet, OpVecBndPressure)>,
    pub oper_fcnt_vel_itf: Vec<OpVecItfContinuity>,
    pub oper_fcnt_pres_itf: Vec<OpSclItfContinuity>,
}

impl TransientMassFlow {
    pub fn new(num_comp: usize) -> TransientMassFlow {
        // initialize struct
        let mut phys = TransientMassFlow::default();

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
        phys.mout_bnd = vec![Vec::new(); num_comp];
        phys.mcnt_itf = vec![Vec::new(); num_comp];
        phys.mcnt_lmd = vec![HashMap::new(); num_comp];
        phys.mres_itf = vec![Vec::new(); num_comp];
        phys.mres_mtrn = vec![HashMap::new(); num_comp];

        // return struct
        phys
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

    pub fn add_mout_bnd(&mut self, comp_i: usize, bnd_id: usize) {
        // boundary domain must also be a flow domain
        self.mout_bnd[comp_i].push(bnd_id);
    }

    pub fn add_mcnt_itf(&mut self, comp_i: usize, itf_id: usize, lmd_id: usize) {
        self.mcnt_itf[comp_i].push(itf_id);
        self.mcnt_lmd[comp_i].insert(itf_id, lmd_id);
    }

    pub fn add_mres_itf(&mut self, comp_i: usize, itf_id: usize, mtrn_id: usize) {
        self.mres_itf[comp_i].push(itf_id);
        self.mres_mtrn[comp_i].insert(itf_id, mtrn_id);
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

impl TransientBase for TransientMassFlow {
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

        // mass unknowns
        for comp_i in 0..self.num_comp {
            for (&dom_id, &conc_id) in &self.itr_conc[comp_i] {
                let num_node = vars.dom[dom_id].num_node;
                vars.scl_dom[conc_id].scl_type = ScalarDomainType::Unknown { start: xid };
                xid += num_node;
            }
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
        let mut vel_dir_nid: HashSet<usize> = HashSet::new();
        let mut pres_dir_nid: HashSet<usize> = HashSet::new();

        // iterate over velocity and pressure boundaries
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

        // iterate over mass scalars and flag dirichlet boundaries
        // flag per domain so a Dirichlet BC on one domain does not mark the
        // shared interface node as Dirichlet in the neighboring domain
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
        // mass internal operators
        for comp_i in 0..self.num_comp {
            for &dom_id in &self.itr_mdom {
                // time operator -> comp_i
                let conc_id = self.itr_conc[comp_i][&dom_id];
                let oper_time = OpSclDomTimeUnity::new(dom_id, conc_id);
                self.oper_itr_time.push(oper_time);

                // source operator -> comp_i
                let msrc_id = self.itr_msrc[comp_i][&dom_id];
                let oper_src = OpSclDomSource::new(dom_id, msrc_id, conc_id);
                self.oper_itr_msrc.push(oper_src);

                // diffusion operator -> (comp_i, comp_j)
                for (&comp_j, &diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let conc_j_id = self.itr_conc[comp_j][&dom_id];
                    let oper_diff = OpSclDomDiffusion::new(dom_id, diff_id, conc_id, conc_j_id);
                    self.oper_itr_diff.push(oper_diff);
                }

                // advective mass transport on overlapping mass+flow domains
                if self.itr_vel.contains_key(&dom_id) {
                    let vel_id = self.itr_vel[&dom_id];
                    let tau_diff_id = self.itr_diff[comp_i][&dom_id][&comp_i];
                    let mut diff_drv_ids = Vec::new();
                    for (&comp_j, &diff_id) in &self.itr_diff[comp_i][&dom_id] {
                        let conc_j_id = self.itr_conc[comp_j][&dom_id];
                        diff_drv_ids.push((diff_id, conc_j_id));
                    }
                    let oper_adv = OpSclDomAdvectionUnity::new(dom_id, vel_id, conc_id);
                    let oper_supg = OpSclDomSupgSteadyUnity::new(dom_id, tau_diff_id, vel_id, msrc_id, conc_id, diff_drv_ids);
                    let oper_supg_time = OpSclDomSupgTimeUnity::new(dom_id, tau_diff_id, vel_id, conc_id);
                    self.oper_itr_mass_adv.push((oper_adv, oper_supg, oper_supg_time));
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

            // boundary outflow operator
            for &bnd_id in &self.mout_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let vel_id = self.itr_vel[&dom_id];
                let oper_out = OpSclBndOutflowUnity::new(bnd_id, vel_id, dom_conc_id);
                self.oper_bnd_mout.push(oper_out);
            }

            // interface continuity operator
            for &itf_id in &self.mcnt_itf[comp_i] {
                let dom1_id = vars.itf[itf_id].dom1_id;
                let dom2_id = vars.itf[itf_id].dom2_id;
                let dom_conc1_id = self.itr_conc[comp_i][&dom1_id];
                let dom_conc2_id = self.itr_conc[comp_i][&dom2_id];
                let lmd_id = self.mcnt_lmd[comp_i][&itf_id];
                let oper_cont = OpSclItfContinuity::new(itf_id, lmd_id, dom_conc1_id, dom_conc2_id);
                self.oper_mcnt_itf.push(oper_cont);
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
            let oper_den_time = OpSclDomDensityTime::new(dom_id, den_id, pres_id);
            let oper_div = OpSclDomDivergence::new(dom_id, den_id, vel_id, pres_id);
            let oper_pspg = OpSclDomPspgSteady::new(dom_id, den_id, visc_id, vel_id, pres_id, fce_id, pres_id);
            let oper_pspg_time = OpSclDomPspgTime::new(dom_id, den_id, visc_id, vel_id, pres_id);

            self.oper_itr_flow.push((oper_time, oper_adv, oper_pres, oper_diff, oper_src, oper_supg, oper_supg_time, oper_den_time, oper_div, oper_pspg, oper_pspg_time));
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
        for (oper_time, oper_adv, oper_pres, oper_diff, oper_src, oper_supg, oper_supg_time, oper_den_time, oper_div, oper_pspg, oper_pspg_time) in &self.oper_itr_flow {
            oper_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
            oper_adv.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_pres.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_diff.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_src.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_supg.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_supg_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
            oper_den_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
            oper_div.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_pspg.apply(vars, a_triplet, b_vec, t, 1.0);
            oper_pspg_time.apply_time(vars, a_triplet, b_vec, t, dt, 1.0);
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
        for (oper_adv, oper_supg, oper_supg_time) in &self.oper_itr_mass_adv {
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
        for oper_out in &self.oper_bnd_mout {
            oper_out.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble flow interface data
        for oper_cont in &self.oper_fcnt_vel_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_cont in &self.oper_fcnt_pres_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }

        // assemble mass interface data
        for oper_cont in &self.oper_mcnt_itf {
            oper_cont.apply(vars, a_triplet, b_vec, t, 1.0);
        }
        for oper_mres in &self.oper_mres_itf {
            oper_mres.apply(vars, a_triplet, b_vec, t, 1.0);
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
