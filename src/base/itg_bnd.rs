use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::shape::prelude::*;

#[derive(Default)]
pub struct IntegralBoundary {
    // ids
    pub itgbnd_id: usize,
    pub bnd_id: usize, // boundary this integral is attached to

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

impl IntegralBoundary {
    pub fn new(itgbnd_id: usize, bnd: &Boundary) -> Result<IntegralBoundary, FEChemError> {
        // create struct
        let mut itgbnd = IntegralBoundary::default();
        itgbnd.itgbnd_id = itgbnd_id;
        itgbnd.bnd_id = bnd.bnd_id;

        // iterate through quadrature points
        for eid in 0..bnd.num_elem {
            // get element-specific data
            let (num_quad, num_node, quad_w, quad_n, quad_gna) = match bnd.elem_type[eid] {
                2 => {
                    let num_quad = Lin2.num_quad();
                    let num_node = Lin2.num_node();
                    let quad_w = Lin2.w();
                    let (quad_n, quad_gna) = Lin2.n();
                    (num_quad, num_node, quad_w, quad_n, quad_gna)
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
            let node_id = &bnd.elem_node_id[eid];
            let mut node_x = vec![0.0; num_node];
            let mut node_y = vec![0.0; num_node];
            for i in 0..num_node {
                let nid = node_id[i];
                node_x[i] = bnd.node_x[nid];
                node_y[i] = bnd.node_y[nid];
            }

            // iterate through quadrature points
            for qid in 0..num_quad {
                // shape functions
                let n = &quad_n[qid];
                let dnda = &quad_gna[qid];

                // physical quadrature coordinates
                for i in 0..num_node {
                    quad_x[qid] += n[i] * node_x[i];
                    quad_y[qid] += n[i] * node_y[i];
                }

                
                // tangent: [dx/da, dy/da]; inward normal: [-dy/da, dx/da]
                // J = [dx/da  -dy/da]
                //     [dy/da   dx/da]
                let mut dxda = 0.0;
                let mut dyda = 0.0;
                for i in 0..num_node {
                    dxda += dnda[i] * node_x[i];
                    dyda += dnda[i] * node_y[i];
                }
                jac_mat[qid] = [[dxda, -dyda], [dyda, dxda]];

                // det(J)
                let det = dxda * dxda + dyda * dyda;
                jac_det[qid] = det;

                // J^{-1} = [da/dx  da/dy]
                //          [db/dx  db/dy]
                let inv = [[dxda / det, dyda / det], [dyda / det, dxda / det]];
                jac_inv[qid] = inv;

                // physical gradients
                // grad_x N = J^{-T} grad_ref N
                // [dN/dx] = [da/dx  db/dx] [dN/da]
                // [dN/dy]   [da/dy  db/dy] [dN/db]
                for i in 0..num_node {
                    quad_gnx[qid][i] = inv[0][0] * dnda[i];
                    quad_gny[qid][i] = inv[0][1] * dnda[i];
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
            itgbnd.num_quad.push(num_quad);
            itgbnd.quad_x.push(quad_x);
            itgbnd.quad_y.push(quad_y);
            itgbnd.quad_w.push(quad_w);
            itgbnd.quad_n.push(quad_n);
            itgbnd.quad_gnx.push(quad_gnx);
            itgbnd.quad_gny.push(quad_gny);
            itgbnd.jac_mat.push(jac_mat);
            itgbnd.jac_inv.push(jac_inv);
            itgbnd.jac_met.push(jac_met);
            itgbnd.jac_det.push(jac_det);

        }

        // result
        Ok(itgbnd)
    }
}
