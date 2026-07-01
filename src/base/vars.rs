use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::itg_bnd::IntegralBoundary;
use crate::base::itg_dom::IntegralDomain;
use crate::base::mesh::Mesh;
use crate::base::scl_bnd::ScalarBoundary;
use crate::base::scl_dom::{ScalarDomain, ScalarDomainType};
use faer::Col;

#[derive(Default)]
pub struct Variables {
    // geometry
    pub mesh: Mesh,
    pub dom: Vec<Domain>,
    pub bnd: Vec<Boundary>,

    // integrals
    pub itg_dom: Vec<IntegralDomain>,
    pub itg_bnd: Vec<IntegralBoundary>,

    // scalars
    pub scl_dom: Vec<ScalarDomain>,
    pub scl_bnd: Vec<ScalarBoundary>,
}

impl Variables {
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

    pub fn add_scldom_con(&mut self, dom_id: usize, value_const: f64, file_path: String) -> Result<usize, FEChemError> {
        let scldom_id = self.scl_dom.len();
        let scldom = ScalarDomain::new_from_constant(scldom_id, &self.dom[dom_id], value_const, file_path)?;
        self.scl_dom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_unk(&mut self, dom_id: usize, value_init: f64, file_path: String) -> Result<usize, FEChemError> {
        let scldom_id = self.scl_dom.len();
        let scldom = ScalarDomain::new_from_unknown(scldom_id, &self.dom[dom_id], value_init, file_path)?;
        self.scl_dom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_sclbnd_con(&mut self, bnd_id: usize, value_const: f64) -> Result<usize, FEChemError> {
        let sclbnd_id = self.scl_bnd.len();
        let sclbnd = ScalarBoundary::new_from_constant(sclbnd_id, &self.bnd[bnd_id], value_const)?;
        self.scl_bnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn update_unknown(&mut self, x_vec: &Col<f64>) {
        // iterate over unknown scalars
        for scldom in self.scl_dom.iter_mut() {
            if scldom.scl_type == ScalarDomainType::Unknown {
                let dom = &self.dom[scldom.dom_id];
                let num_node = dom.num_node;
                for nid in 0..num_node {
                    let xid = scldom.unk_start + nid;
                    let value = x_vec[xid];
                    scldom.node_value[nid] = value;
                }
            }
        }
    }

    pub fn write_scalar(&self, ts: usize) -> Result<(), FEChemError> {
        // iterate over writers
        for scldom in self.scl_dom.iter() {
            let dom = &self.dom[scldom.dom_id];
            scldom.write(dom, ts)?;
        }

        Ok(())
    }

}
