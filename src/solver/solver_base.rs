use faer::Col;
use faer::sparse::SparseColMat;

pub trait SolverBase {
    fn solve(&self, a_mat: &SparseColMat<usize, f64>, b_vec: &Col<f64>, x_init: &Col<f64>, mat_size: usize) -> Col<f64>;
}
