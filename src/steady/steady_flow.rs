use crate::base::scl_dom::ScalarDomainType;
use crate::base::vec_dom::VectorDomainType;
use crate::base::vars::Variables;
use crate::operator::prelude::*;
use crate::steady::steady_base::SteadyBase;
use faer::Col;
use faer::sparse::{SparseColMat, Triplet};
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct SteadyFlow {
    // internal data
    pub itr_dom: Vec<usize>,             // dom id
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

    // TODO: add interface

    // reference pressure
    pub pref_use: bool,
    pub pref_dom: usize,
    pub pref_nid: usize,
    pub pref_pres: f64,

    // operators
    pub oper_itr: Vec<(OpVecDomAdvection, OpVecDomPressure, OpVecDomDiffusion, OpVecDomSource, OpVecDomSupg, OpSclDomDivergence, OpSclDomPspg)>,
    pub oper_bnd_vel: Vec<OpVecBndDirichlet>,
    pub oper_bnd_pres: Vec<OpSclBndDirichlet>,

}

impl SteadyFlow {
    pub fn new() -> SteadyFlow {
        SteadyFlow::default()
    }

    pub fn add_flow_dom(&mut self, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) {
        self.itr_dom.push(dom_id);
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

    pub fn set_pres_ref(&mut self, dom_id: usize, nid: usize, pres: f64) {
        self.pref_use = true;
        self.pref_dom = dom_id;
        self.pref_nid = nid;
        self.pref_pres = pres;
    }
}

impl SteadyBase for SteadyFlow {
    fn assemble_operator(&mut self, vars: &mut Variables, mat_size: &mut usize) {
        // step 1: compute matrix size

        // assign start indices for unknowns
        let mut xid = 0;
        for (&dom_id, &vel_id) in &self.itr_vel {
            let num_node = vars.dom[dom_id].num_node;
            vars.vec_dom[vel_id].vec_type = VectorDomainType::Unknown { start: xid };
            xid += 2*num_node;
        }
        for (&dom_id, &pres_id) in &self.itr_pres {
            let num_node = vars.dom[dom_id].num_node;
            vars.scl_dom[pres_id].scl_type = ScalarDomainType::Unknown { start: xid };
            xid += num_node;
        }

        // set matrix size to last unknown index
        *mat_size = xid;

        // step 2: flag dirichlet boundaries

        // list of mesh node ids with dirichlet boundaries
        let mut vel_dir_nid: HashSet<usize> = HashSet::new();
        let mut pres_dir_nid: HashSet<usize> = HashSet::new();

        // iterate over velocity and pressureboundaries
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

        // iterate over scalars and flag dirichlet boundaries
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
        
        // step 3: assemble operators

        // internal operator
        for &dom_id in &self.itr_dom {
            let vel_id = self.itr_vel[&dom_id];
            let pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let visc_id = self.itr_visc[&dom_id];
            let fce_id = self.itr_fce[&dom_id];

            let oper_adv = OpVecDomAdvection::new(dom_id, den_id, vel_id, vel_id);
            let oper_pres = OpVecDomPressure::new(dom_id, vel_id, pres_id);
            let oper_diff = OpVecDomDiffusion::new(dom_id, visc_id, vel_id, vel_id);
            let oper_src = OpVecDomSource::new(dom_id, fce_id, vel_id);
            let oper_supg = OpVecDomSupg::new(dom_id, den_id, visc_id, vel_id, pres_id, fce_id, vel_id);
            let oper_div = OpSclDomDivergence::new(dom_id, den_id, vel_id, pres_id);
            let oper_pspg = OpSclDomPspg::new(dom_id, den_id, visc_id, vel_id, pres_id, fce_id, pres_id);

            self.oper_itr.push((oper_adv, oper_pres, oper_diff, oper_src, oper_supg, oper_div, oper_pspg));
        }

        // boundary velocity operator
        for &bnd_id in &self.vel_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let bnd_vel_id = self.vel_vel[&bnd_id];
            let oper_dir = OpVecBndDirichlet::new(bnd_id, bnd_vel_id, dom_vel_id);
            self.oper_bnd_vel.push(oper_dir);
        }

        // boundary pressure operator
        for &bnd_id in &self.pres_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_pres_id = self.itr_pres[&dom_id];
            let bnd_pres_id = self.pres_pres[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(bnd_id, bnd_pres_id, dom_pres_id);
            self.oper_bnd_pres.push(oper_dir);
        }

    }

    fn assemble_matrix(&self, vars: &Variables, a_mat: &mut SparseColMat<usize, f64>, b_vec: &mut Col<f64>, mat_size: usize) {
        // initialize triplet for matrix assembly
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        *b_vec = Col::zeros(mat_size);

        // assemble internal data
        for (oper_adv, oper_pres, oper_diff, oper_src, oper_supg, oper_div, oper_pspg) in &self.oper_itr {
            oper_adv.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
            oper_pres.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
            oper_diff.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
            oper_src.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
            oper_supg.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
            oper_div.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
            oper_pspg.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }

        // assemble boundary data
        for oper_dir in &self.oper_bnd_vel {
            oper_dir.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
        }
        for oper_dir in &self.oper_bnd_pres {
            oper_dir.apply(vars, &mut a_triplet, b_vec, 0.0, 1.0);
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

        // create sparse matrix from triplet
        *a_mat = SparseColMat::try_new_from_triplets(mat_size, mat_size, &a_triplet).expect("Failed to create sparse matrix from triplets.");
    }
}
