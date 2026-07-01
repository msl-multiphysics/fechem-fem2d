use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum FEChemError {

    // file - general
    #[error("{caller}: File format is not supported. Needs {type_need}. Got {type_got}.")]
    UnsupportedFileFormat {
        caller: String,
        type_need: String,
        type_got: String,
    },

    #[error("{caller}: Could not write file {file_path}.")]
    FileWriteError {
        caller: String,
        file_path: String,
    },

    #[error("Invalid element type")]
    InvalidElementType,

    #[error("Boundary edge with domain nodes {node0} and {node1} not found on any domain element")]
    BoundaryEdgeNotFound { node0: usize, node1: usize },

    // physics solver - general
    #[error("{caller}: Need at least one iteration. Got max_iter = {max_iter}.")]
    InvalidMaxIter {
        caller: String,
        max_iter: usize,
    },
    #[error("{caller}: Tolerance must be greater than zero. Got tol = {tol}.")]
    InvalidTolerance {
        caller: String,
        tol: f64,
    },
    #[error("{caller}: Damping factor must be in range (0, 1]. Got damp = {damp}.")]
    InvalidDamping {
        caller: String,
        damp: f64,
    },
    #[error("{caller}: Failed to converge in {max_iter} iterations.")]
    FailedConvergence {
        caller: String,
        max_iter: usize,
    },


    // matrix solver - general
    #[error("{caller}: Number of threads must be positive. Got num_thread = {num_thread}.")]
    InvalidThreadCount {
        caller: String,
        num_thread: usize,
    },
    #[error("{caller}: Could not solve matrix equation.")]
    FailedMatrixSolve {
        caller: String,
    },

}