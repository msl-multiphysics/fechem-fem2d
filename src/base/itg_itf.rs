use crate::base::error::FEChemError;
use crate::base::geom_itf::Interface;
use crate::shape::prelude::*;

#[derive(Default)]
pub struct IntegralInterface {
    // ids
    pub itgitf_id: usize,
    pub itf_id: usize, // interface this integral is attached to

    // indexing
    // e - element; q - quadrature point; v - node point

    // for domain 1
    pub quad1_x: Vec<Vec<f64>>,            // [e][q] -> x coordinates
    pub quad1_y: Vec<Vec<f64>>,            // [e][q] -> y coordinates
    pub quad1_gnx: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt x
    pub quad1_gny: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt y
    pub jac1_mat: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> Jacobian matrix
    pub jac1_inv: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> inverse Jacobian matrix

    // for domain 2
    pub quad2_x: Vec<Vec<f64>>,            // [e][q] -> x coordinates
    pub quad2_y: Vec<Vec<f64>>,            // [e][q] -> y coordinates
    pub quad2_gnx: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt x
    pub quad2_gny: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt y
    pub jac2_mat: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> Jacobian matrix
    pub jac2_inv: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> inverse Jacobian matrix

    // quadrature point data
    // e - element; q - quadrature point; v - node point
    pub num_quad: Vec<usize>, // [e] -> number of quadrature points per element
    pub quad_w: Vec<Vec<f64>>, // [e][q] -> quadrature weights
    pub quad_n: Vec<Vec<Vec<f64>>>, // [e][q][v] -> shape functions
    pub jac_met: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> metric tensor
    pub jac_det: Vec<Vec<f64>>, // [e][q] -> Jacobian determinant

}

impl IntegralInterface {
    pub fn new(itgitf_id: usize, itf: &Interface) -> Result<IntegralInterface, FEChemError> {
        // create struct
        let mut itgitf = IntegralInterface::default();
        itgitf.itgitf_id = itgitf_id;
        itgitf.itf_id = itf.itf_id;

        // iterate through quadrature points
        for eid in 0..itf.num_elem {
            // get element-specific data
            let (num_quad, num_node, quad_w, quad_n, quad_gna) = match itf.elem_type[eid] {
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

            // compute for domain 1

            // initialize per element storage
            let mut quad1_x = vec![0.0; num_quad];
            let mut quad1_y = vec![0.0; num_quad];
            let mut quad1_gnx = vec![vec![0.0; num_node]; num_quad];
            let mut quad1_gny = vec![vec![0.0; num_node]; num_quad];
            let mut jac1_mat = vec![[[0.0; 2]; 2]; num_quad];
            let mut jac1_inv = vec![[[0.0; 2]; 2]; num_quad];
            let mut jac_met = vec![[[0.0; 2]; 2]; num_quad]; // permutation invariant
            let mut jac_det = vec![0.0; num_quad]; // permutation invariant

            // get element nodes
            let node_id = &itf.elem_node1_id[eid];
            let mut node_x = vec![0.0; num_node];
            let mut node_y = vec![0.0; num_node];
            for i in 0..num_node {
                let nid = node_id[i];
                node_x[i] = itf.node_x[nid];
                node_y[i] = itf.node_y[nid];
            }

            // iterate through quadrature points
            for qid in 0..num_quad {
                // shape functions
                let n = &quad_n[qid];
                let dnda = &quad_gna[qid];

                // physical quadrature coordinates
                for i in 0..num_node {
                    quad1_x[qid] += n[i] * node_x[i];
                    quad1_y[qid] += n[i] * node_y[i];
                }

                // tangent: [dx/da, dy/da]; normal: [-dy/da, dx/da]
                // J = [dx/da  -dy/da]
                //     [dy/da   dx/da]
                let mut dxda = 0.0;
                let mut dyda = 0.0;
                for i in 0..num_node {
                    dxda += dnda[i] * node_x[i];
                    dyda += dnda[i] * node_y[i];
                }
                jac1_mat[qid] = [[dxda, -dyda], [dyda, dxda]];

                // det(J)
                let det = dxda * dxda + dyda * dyda;
                jac_det[qid] = det;

                // J^{-1} = [da/dx  da/dy]
                //          [db/dx  db/dy]
                let inv = [[dxda / det, dyda / det], [dyda / det, dxda / det]];
                jac1_inv[qid] = inv;

                // physical gradients
                // grad_x N = J^{-T} grad_ref N
                // [dN/dx] = [da/dx  db/dx] [dN/da]
                // [dN/dy]   [da/dy  db/dy] [dN/db]
                for i in 0..num_node {
                    quad1_gnx[qid][i] = inv[0][0] * dnda[i];
                    quad1_gny[qid][i] = inv[0][1] * dnda[i];
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
            itgitf.num_quad.push(num_quad);
            itgitf.quad1_x.push(quad1_x);
            itgitf.quad1_y.push(quad1_y);
            itgitf.quad1_gnx.push(quad1_gnx);
            itgitf.quad1_gny.push(quad1_gny);
            itgitf.jac1_mat.push(jac1_mat);
            itgitf.jac1_inv.push(jac1_inv);

            // compute for domain 2

            // initialize per element storage
            let mut quad2_x = vec![0.0; num_quad];
            let mut quad2_y = vec![0.0; num_quad];
            let mut gradn2_x = vec![vec![0.0; num_node]; num_quad];
            let mut gradn2_y = vec![vec![0.0; num_node]; num_quad];
            let mut jac2_mat = vec![[[0.0; 2]; 2]; num_quad];
            let mut jac2_inv = vec![[[0.0; 2]; 2]; num_quad];

            // get element nodes
            let node_id = &itf.elem_node2_id[eid];
            let mut node_x = vec![0.0; num_node];
            let mut node_y = vec![0.0; num_node];
            for i in 0..num_node {
                let nid = node_id[i];
                node_x[i] = itf.node_x[nid];
                node_y[i] = itf.node_y[nid];
            }

            // iterate through quadrature points
            for qid in 0..num_quad {
                // shape functions
                let n = &quad_n[qid];
                let dnda = &quad_gna[qid];

                // physical quadrature coordinates
                for i in 0..num_node {
                    quad2_x[qid] += n[i] * node_x[i];
                    quad2_y[qid] += n[i] * node_y[i];
                }

                // tangent: [dx/da, dy/da]; normal: [-dy/da, dx/da]
                // J = [dx/da  -dy/da]
                //     [dy/da   dx/da]
                let mut dxda = 0.0;
                let mut dyda = 0.0;
                for i in 0..num_node {
                    dxda += dnda[i] * node_x[i];
                    dyda += dnda[i] * node_y[i];
                }
                jac2_mat[qid] = [[dxda, -dyda], [dyda, dxda]];

                // det(J)
                let det = dxda * dxda + dyda * dyda;

                // J^{-1} = [da/dx  da/dy]
                //          [db/dx  db/dy]
                let inv = [[dxda / det, dyda / det], [dyda / det, dxda / det]];
                jac2_inv[qid] = inv;

                // physical gradients
                // grad_x N = J^{-T} grad_ref N
                // [dN/dx] = [da/dx  db/dx] [dN/da]
                // [dN/dy]   [da/dy  db/dy] [dN/db]
                for i in 0..num_node {
                    gradn2_x[qid][i] = inv[0][0] * dnda[i];
                    gradn2_y[qid][i] = inv[0][1] * dnda[i];
                }
            }

            // store element data
            itgitf.quad2_x.push(quad2_x);
            itgitf.quad2_y.push(quad2_y);
            itgitf.quad2_gnx.push(gradn2_x);
            itgitf.quad2_gny.push(gradn2_y);
            itgitf.jac2_mat.push(jac2_mat);
            itgitf.jac2_inv.push(jac2_inv);
 
            // store permutation invariant data
            itgitf.quad_w.push(quad_w);
            itgitf.quad_n.push(quad_n);
            itgitf.jac_met.push(jac_met);
            itgitf.jac_det.push(jac_det);

        }

        // result
        Ok(itgitf)
    }
}
