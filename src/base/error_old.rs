use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum FVChemError {

    // file - general
    #[error("{caller}: File format is not supported. Needs {type_need}. Got {type_got}.")]
    UnsupportedFileFormat {
        caller: String,
        type_need: String,
        type_got: String,
    },
    #[error("{caller}: Could not write to file: {file_path}.")]
    FileWriteError {
        caller: String,
        file_path: String,
    },

    // geometry - bounds
    #[error("{caller}: Minimum x bound must be less than maximum. Got x_min = {x_min}; x_max = {x_max}.")]
    InvalidBoundsX {
        caller: String,
        x_min: f64,
        x_max: f64,
    },
    #[error("{caller}: Minimum y bound must be less than maximum. Got y_min = {y_min}; y_max = {y_max}.")]
    InvalidBoundsY {
        caller: String,
        y_min: f64,
        y_max: f64,
    },
    #[error("{caller}: Needs at least two elements along x. Got x_num = {x_num}.")]
    InvalidNodeCountX {
        caller: String,
        x_num: usize,
    },
    #[error("{caller}: Needs at least two elements along y. Got y_num = {y_num}.")]
    InvalidNodeCountY {
        caller: String,
        y_num: usize,
    },

    // geometry - mesh
    #[error("{caller}: Failed to read mesh file '{file_path}': {source}.")]
    MeshFileRead {
        caller: String,
        file_path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("{caller}: Invalid GMSH mesh: {message}")]
    InvalidGmsh {
        caller: String,
        message: String,
    },
    #[error("{caller}: Interface face not found on dom1d_b: mesh face {gfid}.")]
    InvalidInterface {
        caller: String,
        gfid: usize,
    },
    #[error("{caller}: Interface dom1d face counts do not match. dom1d_a has {num_face_a}; dom1d_b has {num_face_b}.")]
    InvalidInterfaceCount {
        caller: String,
        num_face_a: usize,
        num_face_b: usize,
    },
    
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

    // physics solver - transient
    #[error("{caller}: Time step size must be greater than zero. Got dt = {dt}.")]
    InvalidTimeStep {
        caller: String,
        dt: f64,
    },
    #[error("{caller}: Need at least one time step. Got num_ts = {num_ts}.")]
    InvalidNumTimeSteps {
        caller: String,
        num_ts: usize,
    },
    #[error("{caller}: Write frequency must be at least one. Got num_ts_write = {num_ts_write}.")]
    InvalidWriteFrequency {
        caller: String,
        num_ts_write: usize,
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
