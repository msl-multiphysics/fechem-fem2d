use crate::base::error::FEChemError;
use crate::base::geom_dom::Domain;
use crate::shape::prelude::*;

#[derive(Default)]
pub struct IntegralDomain {
    // ids
    pub itgdom_id: usize, // must be the same as the domain id
    pub dom_id: usize,    // domain this integral is attached to

    // indexing
    // e - element; q - quadrature point; v - node point

    // quadrature points
    pub num_quad: Vec<usize>, // [e] -> number of quadrature points per element
    pub quad_x: Vec<Vec<f64>>, // [e][q] -> x coordinates
    pub quad_y: Vec<Vec<f64>>, // [e][q] -> y coordinates

    // shape functions
    pub quad_w: Vec<Vec<f64>>, // [e][q] -> quadrature weights
    pub quad_n: Vec<Vec<Vec<f64>>>, // [e][q][v] -> shape functions
    pub quad_gnx: Vec<Vec<Vec<f64>>>, // [e][q][v] -> gradient of shape function wrt x
    pub quad_gny: Vec<Vec<Vec<f64>>>, // [e][q][v] -> gradient of shape function wrt y

    // Jacobian matrix
    pub jac_mat: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> Jacobian matrix
    pub jac_inv: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> inverse Jacobian matrix
    pub jac_met: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> metric tensor
    pub jac_det: Vec<Vec<f64>>, // [e][q] -> Jacobian determinant
}

impl IntegralDomain {
    pub fn new(itgdom_id: usize, dom: &Domain) -> Result<IntegralDomain, FEChemError> {
        // create struct
        let mut itgdom = IntegralDomain::default();
        itgdom.itgdom_id = itgdom_id;
        itgdom.dom_id = dom.dom_id;

        // iterate through elements
        for eid in 0..dom.num_elem {
            // get element-specific data
            let (num_quad, num_node, quad_w, quad_n, quad_gna, quad_gnb) = match dom.elem_type[eid] {
                13 => {
                    let num_quad = Tri3.num_quad();
                    let num_node = Tri3.num_node();
                    let quad_w = Tri3.w();
                    let (quad_n, quad_gna, quad_gnb) = Tri3.n();
                    (num_quad, num_node, quad_w, quad_n, quad_gna, quad_gnb)
                }
                14 => {
                    let num_quad = Quad4.num_quad();
                    let num_node = Quad4.num_node();
                    let quad_w = Quad4.w();
                    let (quad_n, quad_gna, quad_gnb) = Quad4.n();
                    (num_quad, num_node, quad_w, quad_n, quad_gna, quad_gnb)
                }
                _ => {
                    return Err(FEChemError::InvalidElementType);
                }
            };

            // initialize per element storage
            let mut quad_x = vec![0.0; num_quad];
            let mut quad_y = vec![0.0; num_quad];
            let mut quad_gnx = vec![vec![0.0; num_node]; num_quad];
            let mut quad_gny = vec![vec![0.0; num_node]; num_quad];
            let mut jac_mat = vec![[[0.0; 2]; 2]; num_quad];
            let mut jac_inv = vec![[[0.0; 2]; 2]; num_quad];
            let mut jac_met = vec![[[0.0; 2]; 2]; num_quad];
            let mut jac_det = vec![0.0; num_quad];

            // get element nodes
            let node_id = &dom.elem_node_id[eid];
            let mut node_x = vec![0.0; num_node];
            let mut node_y = vec![0.0; num_node];
            for i in 0..num_node {
                let nid = node_id[i];
                node_x[i] = dom.node_x[nid];
                node_y[i] = dom.node_y[nid];
            }

            // iterate through quadrature points
            for qid in 0..num_quad {
                // shape functions
                let n = &quad_n[qid];
                let dnda = &quad_gna[qid];
                let dndb = &quad_gnb[qid];

                // quadrature point coordinates
                for i in 0..num_node {
                    quad_x[qid] += n[i] * node_x[i];
                    quad_y[qid] += n[i] * node_y[i];
                }

                // J = [dx/da  dx/db]
                //     [dy/da  dy/db]
                let mut dxda = 0.0;
                let mut dxdb = 0.0;
                let mut dyda = 0.0;
                let mut dydb = 0.0;
                for i in 0..num_node {
                    dxda += dnda[i] * node_x[i];
                    dxdb += dndb[i] * node_x[i];
                    dyda += dnda[i] * node_y[i];
                    dydb += dndb[i] * node_y[i];
                }
                jac_mat[qid] = [[dxda, dxdb], [dyda, dydb]];

                // det(J)
                let det = dxda * dydb - dxdb * dyda;
                jac_det[qid] = det;

                // J^{-1} = [da/dx  da/dy]
                //          [db/dx  db/dy]
                let inv = [[dydb / det, -dxdb / det], [-dyda / det, dxda / det]];

                jac_inv[qid] = inv;

                // physical gradients
                // grad_x N = J^{-T} grad_ref N
                // [dN/dx] = [da/dx  db/dx] [dN/da]
                // [dN/dy]   [da/dy  db/dy] [dN/db]
                for i in 0..num_node {
                    quad_gnx[qid][i] = inv[0][0] * dnda[i] + inv[1][0] * dndb[i];
                    quad_gny[qid][i] = inv[0][1] * dnda[i] + inv[1][1] * dndb[i];
                }

                // metric tensor
                // G = J^{-1} J^{-T}
                // G_ij = sum_k inv[i][k] * inv[j][k]
                let mut g = [[0.0; 2]; 2];
                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..2 {
                            g[i][j] += inv[i][k] * inv[j][k];
                        }
                    }
                }
                jac_met[qid] = g;
            }

            // store element data
            itgdom.num_quad.push(num_quad);
            itgdom.quad_x.push(quad_x);
            itgdom.quad_y.push(quad_y);
            itgdom.quad_w.push(quad_w);
            itgdom.quad_n.push(quad_n);
            itgdom.quad_gnx.push(quad_gnx);
            itgdom.quad_gny.push(quad_gny);
            itgdom.jac_mat.push(jac_mat);
            itgdom.jac_inv.push(jac_inv);
            itgdom.jac_met.push(jac_met);
            itgdom.jac_det.push(jac_det);

        }

        // result
        Ok(itgdom)
    }
}
