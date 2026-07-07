use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::geom_itf::Interface;
use crate::base::itg_bnd::IntegralBoundary;
use crate::base::itg_dom::IntegralDomain;
use crate::base::itg_itf::IntegralInterface;
use crate::base::mesh::Mesh;
use crate::base::scl_bnd::{ScalarBoundary, update_function_sclbnd};
use crate::base::scl_dom::{ScalarDomain, update_function_scldom};
use crate::base::scl_itf::ScalarInterface;
use crate::base::vec_bnd::{VectorBoundary, update_function_vecbnd};
use crate::base::vec_dom::{VectorDomain, update_function_vecdom};
use crate::base::vec_itf::VectorInterface;
use faer::Col;

#[derive(Default)]
pub struct Variables {
    // geometry
    pub mesh: Mesh,
    pub dom: Vec<Domain>,
    pub bnd: Vec<Boundary>,
    pub itf: Vec<Interface>,

    // integrals
    pub itg_dom: Vec<IntegralDomain>,
    pub itg_bnd: Vec<IntegralBoundary>,
    pub itg_itf: Vec<IntegralInterface>,

    // scalars
    pub scl_dom: Vec<ScalarDomain>,
    pub scl_bnd: Vec<ScalarBoundary>,
    pub scl_itf: Vec<ScalarInterface>,

    // vectors
    pub vec_dom: Vec<VectorDomain>,
    pub vec_bnd: Vec<VectorBoundary>,
    pub vec_itf: Vec<VectorInterface>,
}

impl Variables {
    pub fn new(file_path: String) -> Result<Variables, FEChemError> {
        let mut vars = Variables::default();
        vars.mesh = Mesh::new(file_path)?;
        Ok(vars)
    }

    pub fn new_from_bounds(x_min: f64, y_min: f64, x_max: f64, y_max: f64, num_elem_x: usize, num_elem_y: usize) -> Result<Variables, FEChemError> {
        let mut vars = Variables::default();
        vars.mesh = Mesh::new_from_bounds(x_min, y_min, x_max, y_max, num_elem_x, num_elem_y)?;
        Ok(vars)
    }

    pub fn add_dom(&mut self, reg_id: usize) -> Result<usize, FEChemError> {
        let dom_id = self.dom.len();
        let itgdom_id = self.itg_dom.len();
        let dom = Domain::new(dom_id, &self.mesh, reg_id)?;
        let itgdom = IntegralDomain::new(itgdom_id, &dom)?;
        self.dom.push(dom);
        self.itg_dom.push(itgdom);
        Ok(dom_id)
    }

    pub fn add_bnd(&mut self, dom_id: usize, reg_id: usize) -> Result<usize, FEChemError> {
        let bnd_id = self.bnd.len();
        let itgbnd_id = self.itg_bnd.len();
        let bnd = Boundary::new(bnd_id, &self.mesh, &self.dom[dom_id], reg_id)?;
        let itgbnd = IntegralBoundary::new(itgbnd_id, &bnd)?;
        self.bnd.push(bnd);
        self.itg_bnd.push(itgbnd);
        Ok(bnd_id)
    }

    pub fn add_itf(&mut self, dom1_id: usize, dom2_id: usize, reg_id: usize) -> Result<usize, FEChemError> {
        let itf_id = self.itf.len();
        let itgitf_id = self.itg_itf.len();
        let itf = Interface::new(itf_id, &self.mesh, &self.dom[dom1_id], &self.dom[dom2_id], reg_id)?;
        let itgitf = IntegralInterface::new(itgitf_id, &itf)?;
        self.itf.push(itf);
        self.itg_itf.push(itgitf);
        Ok(itf_id)
    }

    pub fn add_scldom_con(&mut self, dom_id: usize, value_const: f64, file_path: String) -> Result<usize, FEChemError> {
        let scldom_id = self.scl_dom.len();
        let scldom = ScalarDomain::new_from_constant(scldom_id, &self.dom[dom_id], value_const, file_path)?;
        self.scl_dom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_fun(&mut self, dom_id: usize, value_func: Box<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<usize, FEChemError> {
        // TODO: check that scldom_ids are unknown-type scalars
        // function works only on unknown scalars
        let scldom_id = self.scl_dom.len();
        let scldom = ScalarDomain::new_from_function(scldom_id, &self.dom[dom_id], value_func, scldom_ids, file_path)?;
        self.scl_dom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_unk(&mut self, dom_id: usize, value_init: f64, file_path: String) -> Result<usize, FEChemError> {
        let scldom_id = self.scl_dom.len();
        let scldom = ScalarDomain::new_from_unknown(scldom_id, &self.dom[dom_id], value_init, file_path)?;
        self.scl_dom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_sclbnd_con(&mut self, bnd_id: usize, value_const: f64, file_path: String) -> Result<usize, FEChemError> {
        let sclbnd_id = self.scl_bnd.len();
        let sclbnd = ScalarBoundary::new_from_constant(sclbnd_id, &self.bnd[bnd_id], value_const, file_path)?;
        self.scl_bnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclbnd_fun(&mut self, bnd_id: usize, value_func: Box<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<usize, FEChemError> {
        // TODO: check that scldom_ids are unknown-type scalars
        let sclbnd_id = self.scl_bnd.len();
        let sclbnd = ScalarBoundary::new_from_function(sclbnd_id, &self.bnd[bnd_id], value_func, scldom_ids, file_path)?;
        self.scl_bnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclitf_unk(&mut self, itf_id: usize, value_init: f64, file_path: String) -> Result<usize, FEChemError> {
        let sclitf_id = self.scl_itf.len();
        let sclitf = ScalarInterface::new_from_unknown(sclitf_id, &self.itf[itf_id], value_init, file_path)?;
        self.scl_itf.push(sclitf);
        Ok(sclitf_id)
    }

    pub fn add_vecdom_con(&mut self, dom_id: usize, value_const_x: f64, value_const_y: f64, file_path: String) -> Result<usize, FEChemError> {
        let vecdom_id = self.vec_dom.len();
        let vecdom = VectorDomain::new_from_constant(vecdom_id, &self.dom[dom_id], value_const_x, value_const_y, file_path)?;
        self.vec_dom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecdom_fun(&mut self, dom_id: usize, value_func: Box<dyn Fn(f64, &[f64]) -> (f64, f64) + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<usize, FEChemError> {
        let vecdom_id = self.vec_dom.len();
        let vecdom = VectorDomain::new_from_function(vecdom_id, &self.dom[dom_id], value_func, scldom_ids, file_path)?;
        self.vec_dom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecdom_unk(&mut self, dom_id: usize, value_init_x: f64, value_init_y: f64, file_path: String) -> Result<usize, FEChemError> {
        let vecdom_id = self.vec_dom.len();
        let vecdom = VectorDomain::new_from_unknown(vecdom_id, &self.dom[dom_id], value_init_x, value_init_y, file_path)?;
        self.vec_dom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecbnd_con(&mut self, bnd_id: usize, value_const_x: f64, value_const_y: f64, file_path: String) -> Result<usize, FEChemError> {
        let vecbnd_id = self.vec_bnd.len();
        let vecbnd = VectorBoundary::new_from_constant(vecbnd_id, &self.bnd[bnd_id], value_const_x, value_const_y, file_path)?;
        self.vec_bnd.push(vecbnd);
        Ok(vecbnd_id)
    }

    pub fn add_vecbnd_fun(&mut self, bnd_id: usize, value_func: Box<dyn Fn(f64, &[f64]) -> (f64, f64) + Send + Sync>, scldom_ids: Vec<usize>, file_path: String) -> Result<usize, FEChemError> {
        let vecbnd_id = self.vec_bnd.len();
        let vecbnd = VectorBoundary::new_from_function(vecbnd_id, &self.bnd[bnd_id], value_func, scldom_ids, file_path)?;
        self.vec_bnd.push(vecbnd);
        Ok(vecbnd_id)
    }

    pub fn update_unknown(&mut self, x_vec: &Col<f64>) {
        // skips if not unknown type
        for scldom in self.scl_dom.iter_mut() {
            let dom = &self.dom[scldom.dom_id];
            scldom.update_unknown(dom, x_vec);
        }
        for sclitf in self.scl_itf.iter_mut() {
            let itf = &self.itf[sclitf.itf_id];
            sclitf.update_unknown(itf, x_vec);
        }
        for vecdom in self.vec_dom.iter_mut() {
            let dom = &self.dom[vecdom.dom_id];
            vecdom.update_unknown(dom, x_vec);
        }
        for vecitf in self.vec_itf.iter_mut() {
            let itf = &self.itf[vecitf.itf_id];
            vecitf.update_unknown(itf, x_vec);
        }
    }

    pub fn update_function(&mut self, t: f64) {
        for scldom_id in 0..self.scl_dom.len() {
            update_function_scldom(self, scldom_id, t);
        }
        for sclbnd_id in 0..self.scl_bnd.len() {
            update_function_sclbnd(self, sclbnd_id, t);
        }
        for vecdom_id in 0..self.vec_dom.len() {
            update_function_vecdom(self, vecdom_id, t);
        }
        for vecbnd_id in 0..self.vec_bnd.len() {
            update_function_vecbnd(self, vecbnd_id, t);
        }
    }

    pub fn update_prev(&mut self) {
        for scldom in self.scl_dom.iter_mut() {
            scldom.update_prev();
        }
        for sclitf in self.scl_itf.iter_mut() {
            sclitf.update_prev();
        }
        for vecdom in self.vec_dom.iter_mut() {
            vecdom.update_prev();
        }
        for vecitf in self.vec_itf.iter_mut() {
            vecitf.update_prev();
        }
    }

    pub fn write_scalar(&mut self, t: f64, ts: usize) -> Result<(), FEChemError> {
        // compute function values at nodes
        self.update_function(t);

        // iterate over writers
        for scldom in self.scl_dom.iter_mut() {
            let dom = &self.dom[scldom.dom_id];
            scldom.write(dom, ts)?;
        }
        for sclbnd in self.scl_bnd.iter_mut() {
            let bnd = &self.bnd[sclbnd.bnd_id];
            sclbnd.write(bnd, ts)?;
        }
        for vecdom in self.vec_dom.iter_mut() {
            let dom = &self.dom[vecdom.dom_id];
            vecdom.write(dom, ts)?;
        }
        for vecbnd in self.vec_bnd.iter_mut() {
            let bnd = &self.bnd[vecbnd.bnd_id];
            vecbnd.write(bnd, ts)?;
        }

        Ok(())
    }
}
