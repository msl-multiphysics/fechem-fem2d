use crate::base::error::FEChemError;
use crate::solver::solver_base::SolverBase;
use faer::sparse::SparseColMat;
use faer::{Col, Mat, Par};
use faer_gmres::gmres;
use std::num::NonZeroUsize;

#[derive(Default)]
pub struct SolverGmres {}

impl SolverBase for SolverGmres {
    fn solve(
        &self,
        a_mat: &SparseColMat<usize, f64>,
        b_vec: &Col<f64>,
        x_init: &Col<f64>,
        mat_size: usize,
    ) -> Result<Col<f64>, FEChemError> {
        let mut b_mat = Mat::<f64>::zeros(b_vec.nrows(), 1);
        let mut x_mat = Mat::<f64>::zeros(x_init.nrows(), 1);
        for i in 0..b_vec.nrows() {
            b_mat[(i, 0)] = b_vec[i];
        }
        for i in 0..x_init.nrows() {
            x_mat[(i, 0)] = x_init[i];
        }

        let max_iter = mat_size.max(1);
        let tol = 1e-8;

        gmres(
            a_mat.as_ref(),
            b_mat.as_ref(),
            x_mat.as_mut(),
            max_iter,
            tol,
            None,
        )
        .map_err(|_| FEChemError::FailedMatrixSolve {
            caller: "SolverGmres::solve".to_string(),
        })?;

        let mut x_vec = Col::<f64>::zeros(x_mat.nrows());
        for i in 0..x_mat.nrows() {
            x_vec[i] = x_mat[(i, 0)];
        }

        Ok(x_vec)
    }
}

impl SolverGmres {
    pub fn new(num_thread: usize) -> Result<Self, FEChemError> {
        // error handling
        if num_thread == 0 {
            return Err(FEChemError::InvalidThreadCount {
                caller: "SolverGmres::new".to_string(),
                num_thread,
            });
        }

        // set number of threads
        faer::set_global_parallelism(Par::Rayon(
            NonZeroUsize::new(num_thread).expect("Number of threads must be positive"),
        ));

        // return
        Ok(SolverGmres::default())
    }
}
