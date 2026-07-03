use crate::base::geom_bnd::Boundary;
use crate::base::error::FEChemError;
use crate::shape::prelude::*;

#[derive(Default)]
pub struct IntegralBoundary {
    // ids
    pub itgbnd_id: usize,
    pub bnd_id: usize,  // boundary this integral is attached to

    // quadrature point data
    // e - element; q - quadrature point; v - node point
    pub num_quad: Vec<usize>,  // [e] -> number of quadrature points per element
    pub quad_x: Vec<Vec<f64>>,  // [e][q] -> x coordinates
    pub quad_y: Vec<Vec<f64>>,  // [e][q] -> y coordinates
    pub gradn_x: Vec<Vec<Vec<f64>>>,  // [e][q][v] -> gradient of shape function wrt x
    pub gradn_y: Vec<Vec<Vec<f64>>>,  // [e][q][v] -> gradient of shape function wrt y
    pub jac_mat: Vec<Vec<[[f64; 2]; 2]>>,  // [e][q][i][j] -> Jacobian matrix
    pub jac_inv: Vec<Vec<[[f64; 2]; 2]>>,  // [e][q][i][j] -> inverse Jacobian matrix
    pub jac_met: Vec<Vec<[[f64; 2]; 2]>>,  // [e][q][i][j] -> metric tensor
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
            match bnd.elem_node_num[eid] {
                2 => {compute_lin2(&mut itgbnd, bnd, eid);}
                _ => {return Err(FEChemError::InvalidElementType);}
            }
        }

        // result
        Ok(itgbnd)
    }
}

fn compute_lin2(itgbnd: &mut IntegralBoundary, bnd: &Boundary, eid: usize) {
    // number of quadrature points and nodes
    let num_quad = W_LIN2.len();
    let num_node = 2;

    // initialize per element storage
    let mut quad_x = vec![0.0; num_quad];
    let mut quad_y = vec![0.0; num_quad];
    let mut gradn_x = vec![vec![0.0; num_node]; num_quad];
    let mut gradn_y = vec![vec![0.0; num_node]; num_quad];
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
        // reference coordinates of quadrature point
        let a = A_LIN2[qid];

        // shape functions and reference gradients
        let n = lin2_eval(a);
        let dnda = lin2_grad(a);

        // physical quadrature coordinates
        for i in 0..num_node {
            quad_x[qid] += n[i] * node_x[i];
            quad_y[qid] += n[i] * node_y[i];
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
        jac_mat[qid] = [
            [ dxda, -dyda],
            [ dyda,  dxda],
        ];

        // jacobian determinant
        let det = dxda * dxda + dyda * dyda;
        jac_det[qid] = det;

        // J^{-1} = [
        //   da/dx  da/dy
        //   db/dx  db/dy
        // ]
        let inv = [
            [dxda / det, dyda / det],
            [dyda / det, dxda / det],
        ];
        jac_inv[qid] = inv;

        // physical gradients
        // grad_x N = J^{-T} grad_ref N
        // [dN/dx] = [da/dx  db/dx] [dN/da]
        // [dN/dy]   [da/dy  db/dy] [dN/db]
        for i in 0..num_node {
            gradn_x[qid][i] = inv[0][0] * dnda[i];
            gradn_y[qid][i] = inv[0][1] * dnda[i];
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
    itgbnd.num_quad.push(num_quad);
    itgbnd.quad_x.push(quad_x);
    itgbnd.quad_y.push(quad_y);
    itgbnd.gradn_x.push(gradn_x);
    itgbnd.gradn_y.push(gradn_y);
    itgbnd.jac_mat.push(jac_mat);
    itgbnd.jac_inv.push(jac_inv);
    itgbnd.jac_met.push(jac_met);
    itgbnd.jac_det.push(jac_det);
}
