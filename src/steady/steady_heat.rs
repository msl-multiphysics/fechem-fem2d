use crate::base::vars::Variables;
use crate::operator::prelude::*;
use crate::steady::steady_base::SteadyBase;
use faer::Col;
use faer::sparse::{SparseColMat, Triplet};
use std::collections::HashMap;

#[derive(Default)]
pub struct SteadyHeat {
    // internal data
    pub itr_dom: Vec<usize>,  // dom id
    pub itr_temp: HashMap<usize, usize>,  // temperature (unknown)
    pub itr_cond: HashMap<usize, usize>,  // thermal conductivity
    pub itr_hsrc: HashMap<usize, usize>,  // heat source

    // boundary data
    pub temp_bnd: Vec<usize>,  // bnd with temperature BC
    pub temp_temp: HashMap<usize, usize>,  // temperature
    pub hflx_bnd: Vec<usize>,  // bnd with heat flux BC
    pub hflx_hflx: HashMap<usize, usize>,  // heat flux

    // operators
    pub oper_itr: Vec<(OperatorDiffusion, OperatorSource)>,
    pub oper_bnd_temp: Vec<OperatorDirichlet>,
    pub oper_bnd_hflx: Vec<OperatorNeumannDiffusion>,
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
    
}

impl SteadyBase for SteadyHeat {
    fn assemble_operator(&mut self, vars: &mut Variables, mat_size: &mut usize) {
        // step 1: compute matrix size

        // assign start indices for unknowns
        let mut xid = 0;
        for (&dom_id, &temp_id) in &self.itr_temp {
            let num_node = vars.dom[dom_id].num_node;
            vars.scl_dom[temp_id].unk_start = xid;
            xid += num_node;
        }

        // set matrix size to last unknown index
        *mat_size = xid;
        
        // step 2: flag dirichlet boundaries

        // iterate over temperature boundaries
        for &bnd_id in &self.temp_bnd {
            // get domain and temperature ids
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];

            // flag dirichlet boundaries
            for &nid in &vars.bnd[bnd_id].node_bnd_dom_id {
                vars.scl_dom[dom_temp_id].node_dir[nid] = true;
            }
        }

        // step 3: assemble operators

        // internal operator
        for &dom_id in &self.itr_dom {
            let temp_id = self.itr_temp[&dom_id];
            let cond_id = self.itr_cond[&dom_id];
            let hsrc_id = self.itr_hsrc[&dom_id];
            let oper_cond = OperatorDiffusion::new(dom_id, cond_id, temp_id, temp_id);
            let oper_src = OperatorSource::new(dom_id, hsrc_id, temp_id);
            self.oper_itr.push((oper_cond, oper_src));
        }

        // boundary temperature operator
        for &bnd_id in &self.temp_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_temp_id = self.temp_temp[&bnd_id];
            let oper_dir = OperatorDirichlet::new(bnd_id, bnd_temp_id, dom_temp_id);
            self.oper_bnd_temp.push(oper_dir);
        }

        // boundary flux operator
        for &bnd_id in &self.hflx_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_hflx_id = self.hflx_hflx[&bnd_id];
            let oper_neu = OperatorNeumannDiffusion::new(bnd_id, bnd_hflx_id, dom_temp_id);
            self.oper_bnd_hflx.push(oper_neu);
        }

    }

    fn assemble_matrix(&self, vars: &Variables, a_mat: &mut SparseColMat<usize, f64>, b_vec: &mut Col<f64>, mat_size: usize) {
        // initialize triplet for matrix assembly
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        *b_vec = Col::zeros(mat_size);
        
        // assemble internal data
        for (oper_cond, oper_src) in &self.oper_itr {
            oper_cond.apply(vars, &mut a_triplet, b_vec, 1.0);
            oper_src.apply(vars, &mut a_triplet, b_vec, 1.0);
        }

        // assemble boundary data
        for oper_dir in &self.oper_bnd_temp {
            oper_dir.apply(vars, &mut a_triplet, b_vec, 1.0);
        }
        for oper_neu in &self.oper_bnd_hflx {
            oper_neu.apply(vars, &mut a_triplet, b_vec, 1.0);
        }

        // create sparse matrix from triplet
        *a_mat = SparseColMat::try_new_from_triplets(mat_size, mat_size, &a_triplet).expect("Failed to create sparse matrix from triplets.");

    }

}
