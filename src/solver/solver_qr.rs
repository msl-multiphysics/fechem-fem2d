use crate::base::error::FEChemError;
use crate::solver::solver_base::SolverBase;
use faer::{Col, Par};
use faer::sparse::SparseColMat;
use faer::sparse::linalg::solvers::{Qr, SymbolicQr};
use faer::linalg::solvers::Solve;
use std::num::NonZeroUsize;

#[derive(Default)]
pub struct SolverQr {}

impl SolverBase for SolverQr {
    fn solve(&self, a_mat: &SparseColMat<usize, f64>, b_vec: &Col<f64>, _: &Col<f64>, _: usize) -> Result<Col<f64>, FEChemError> {
        let symbolic = SymbolicQr::try_new(a_mat.symbolic()).map_err(|_| FEChemError::FailedMatrixSolve {caller: "SolverQr::solve".to_string()})?;
        let qr = Qr::try_new_with_symbolic(symbolic, a_mat.as_ref()).map_err(|_| FEChemError::FailedMatrixSolve {caller: "SolverQr::solve".to_string()})?;
        Ok(qr.solve(b_vec))
    }
}

impl SolverQr {
    pub fn new(num_thread: usize) -> Result<Self, FEChemError> {
        // error handling
        if num_thread == 0 {
            return Err(FEChemError::InvalidThreadCount {caller: "SolverQr::new".to_string(), num_thread: num_thread});
        }

        // set number of threads
        faer::set_global_parallelism(Par::Rayon(NonZeroUsize::new(num_thread).expect("Number of threads must be positive")));
        
        Ok(SolverQr::default())
    }
}
