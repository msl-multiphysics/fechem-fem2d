use crate::base::error::FEChemError;
use crate::base::geom_itf::Interface;
use crate::shape::prelude::*;

#[derive(Default)]
pub struct IntegralInterface {
    // ids
    pub itgitf_id: usize,
    pub itf_id: usize, // interface this integral is attached to

    // quadrature point data
    // e - element; q - quadrature point; v - node point
    pub num_quad: Vec<usize>, // [e] -> number of quadrature points per element
    pub jac_met: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> metric tensor
    pub jac_det: Vec<Vec<f64>>, // [e][q] -> Jacobian determinant

    // for domain 1
    pub quad1_x: Vec<Vec<f64>>,            // [e][q] -> x coordinates
    pub quad1_y: Vec<Vec<f64>>,            // [e][q] -> y coordinates
    pub gradn1_x: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt x
    pub gradn1_y: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt y
    pub jac1_mat: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> Jacobian matrix
    pub jac1_inv: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> inverse Jacobian matrix

    // for domain 2
    pub quad2_x: Vec<Vec<f64>>,            // [e][q] -> x coordinates
    pub quad2_y: Vec<Vec<f64>>,            // [e][q] -> y coordinates
    pub gradn2_x: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt x
    pub gradn2_y: Vec<Vec<Vec<f64>>>,      // [e][q][v] -> gradient of shape function wrt y
    pub jac2_mat: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> Jacobian matrix
    pub jac2_inv: Vec<Vec<[[f64; 2]; 2]>>, // [e][q][i][j] -> inverse Jacobian matrix
}

impl IntegralInterface {
    pub fn new(itgitf_id: usize, itf: &Interface) -> Result<IntegralInterface, FEChemError> {
        // create struct
        let mut itgitf = IntegralInterface::default();
        itgitf.itgitf_id = itgitf_id;
        itgitf.itf_id = itf.itf_id;

        // iterate through quadrature points
        for eid in 0..itf.num_elem {
            match itf.elem_node_num[eid] {
                2 => {
                    compute_lin2(&mut itgitf, itf, eid);
                }
                _ => {
                    return Err(FEChemError::InvalidElementType);
                }
            }
        }

        // result
        Ok(itgitf)
    }
}

fn compute_lin2(itgitf: &mut IntegralInterface, itf: &Interface, eid: usize) {
    // number of quadrature points and nodes
    let num_quad = W_LIN2.len();
    let num_node = 2;

    // compute for domain 1

    // initialize per element storage
    let mut quad1_x = vec![0.0; num_quad];
    let mut quad1_y = vec![0.0; num_quad];
    let mut gradn1_x = vec![vec![0.0; num_node]; num_quad];
    let mut gradn1_y = vec![vec![0.0; num_node]; num_quad];
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
        // reference coordinates of quadrature point
        let a = A_LIN2[qid];

        // shape functions and reference gradients
        let n = lin2_eval(a);
        let dnda = lin2_grad(a);

        // physical quadrature coordinates
        for i in 0..num_node {
            quad1_x[qid] += n[i] * node_x[i];
            quad1_y[qid] += n[i] * node_y[i];
        }

        // Jacobian entries
        // tangent:  [dx/da, dy/da]
        // normal:   [-dy/da, dx/da]  (completes a 2x2 basis for the line in the plane)
        let mut dxda = 0.0;
        let mut dyda = 0.0;
        for i in 0..num_node {
            dxda += dnda[i] * node_x[i];
            dyda += dnda[i] * node_y[i];
        }
        jac1_mat[qid] = [[dxda, -dyda], [dyda, dxda]];

        // jacobian determinant
        let det = dxda * dxda + dyda * dyda;
        jac_det[qid] = det;

        // J^{-1} = [
        //   da/dx  da/dy
        //   db/dx  db/dy
        // ]
        let inv = [[dxda / det, dyda / det], [dyda / det, dxda / det]];
        jac1_inv[qid] = inv;

        // physical gradients
        // grad_x N = J^{-T} grad_ref N
        // [dN/dx] = [da/dx  db/dx] [dN/da]
        // [dN/dy]   [da/dy  db/dy] [dN/db]
        for i in 0..num_node {
            gradn1_x[qid][i] = inv[0][0] * dnda[i];
            gradn1_y[qid][i] = inv[0][1] * dnda[i];
        }

        // metric tensor for GLS
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
    itgitf.gradn1_x.push(gradn1_x);
    itgitf.gradn1_y.push(gradn1_y);
    itgitf.jac1_mat.push(jac1_mat);
    itgitf.jac1_inv.push(jac1_inv);
    itgitf.jac_met.push(jac_met); // permutation invariant
    itgitf.jac_det.push(jac_det); // permutation invariant

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
        // reference coordinates of quadrature point
        let a = A_LIN2[qid];

        // shape functions and reference gradients
        let n = lin2_eval(a);
        let dnda = lin2_grad(a);

        // physical quadrature coordinates
        for i in 0..num_node {
            quad2_x[qid] += n[i] * node_x[i];
            quad2_y[qid] += n[i] * node_y[i];
        }

        // Jacobian entries
        // tangent:  [dx/da, dy/da]
        // normal:   [-dy/da, dx/da]  (completes a 2x2 basis for the line in the plane)
        let mut dxda = 0.0;
        let mut dyda = 0.0;
        for i in 0..num_node {
            dxda += dnda[i] * node_x[i];
            dyda += dnda[i] * node_y[i];
        }
        jac2_mat[qid] = [[dxda, -dyda], [dyda, dxda]];

        // jacobian determinant
        let det = dxda * dxda + dyda * dyda;

        // J^{-1} = [
        //   da/dx  da/dy
        //   db/dx  db/dy
        // ]
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
    itgitf.gradn2_x.push(gradn2_x);
    itgitf.gradn2_y.push(gradn2_y);
    itgitf.jac2_mat.push(jac2_mat);
    itgitf.jac2_inv.push(jac2_inv);
}
