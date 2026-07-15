use crate::base::error::FEChemError;
use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vars::Variables;
use crate::base::vec_dom::VectorDomainType;
use crate::base::vec_itf::VectorInterfaceType;
use crate::solver::solver_base::SolverBase;
use faer::Col;
use faer::sparse::{Pair, SparseColMat, SymbolicSparseColMat, Triplet};
use std::time::{Duration, Instant};


// base trait for steady-state solvers.
pub trait SteadyBase {
    // to be implemented in specific solvers
    fn initial_matrix(&self, vars: &mut Variables) -> usize; // computes matrix size
    fn initial_dirichlet(&self, vars: &mut Variables); // flags dirichlet boundaries
    fn initial_operator(&mut self, vars: &mut Variables); // initializes operators
    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>);

    fn solve(&mut self, vars: &mut Variables, solver: Box<dyn SolverBase>, max_iter: usize, tol: f64, damp: f64) -> Result<(), FEChemError> {
        let time_start = Instant::now();
        println!("Starting steady-state solver.");

        // error handling
        if max_iter == 0 {
            return Err(FEChemError::InvalidMaxIter {caller: "SteadyBase::solve".to_string(), max_iter});
        }
        if tol <= 0.0 {
            return Err(FEChemError::InvalidTolerance {caller: "SteadyBase::solve".to_string(), tol});
        }
        if damp <= 0.0 || damp > 1.0 {
            return Err(FEChemError::InvalidDamping {caller: "SteadyBase::solve".to_string(), damp});
        }

        // initialize time measurement
        let mut time_assemble = Duration::ZERO;
        let mut time_solve = Duration::ZERO;

        let time_0 = Instant::now();

        // initialize operators
        let mat_size = self.initial_matrix(vars);
        self.initial_dirichlet(vars);
        self.initial_operator(vars);

        // initialize solution vectors
        let mut x_udmp_vec: Col<f64>;
        let mut x_iter_vec: Col<f64> = Col::zeros(mat_size);
        
        // load initial unknown values into solution vector
        for scldom in &vars.scl_dom {
            if let ScalarDomainType::Unknown { start } = scldom.scl_type {
                let dom = &vars.dom[scldom.dom_id];
                for nid in 0..dom.num_node {
                    x_iter_vec[start + nid] = scldom.node_value[nid];
                }
            }
        }
        for sclitf in &vars.scl_itf {
            if let ScalarInterfaceType::Unknown { start } = sclitf.scl_type {
                let itf = &vars.itf[sclitf.itf_id];
                for nid in 0..itf.num_node {
                    x_iter_vec[start + nid] = sclitf.node_value[nid];
                }
            }
        }
        for vecdom in &vars.vec_dom {
            if let VectorDomainType::Unknown { start } = vecdom.vec_type {
                let dom = &vars.dom[vecdom.dom_id];
                let num_node = dom.num_node;
                for nid in 0..num_node {
                    x_iter_vec[start + nid] = vecdom.node_value_x[nid];
                    x_iter_vec[start + nid + num_node] = vecdom.node_value_y[nid];
                }
            }
        }
        for vecitf in &vars.vec_itf {
            if let VectorInterfaceType::Unknown { start } = vecitf.vec_type {
                let itf = &vars.itf[vecitf.itf_id];
                let num_node = itf.num_node;
                for nid in 0..num_node {
                    x_iter_vec[start + nid] = vecitf.node_value_x[nid];
                    x_iter_vec[start + nid + num_node] = vecitf.node_value_y[nid];
                }
            }
        }

        // initialize A matrix (triplet form) and b vector
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let mut b_vec: Col<f64> = Col::zeros(mat_size);
        
        // make initial assembly of A and b
        // this will be used in the iteration and for finding the sparsity pattern
        self.assemble_matrix(vars, &mut a_triplet, &mut b_vec);
        let a_pair: Vec<Pair<usize, usize>> = a_triplet.iter().map(|t| Pair::new(t.row, t.col)).collect();
        let (symbolic, argsort) = SymbolicSparseColMat::try_new_from_indices(mat_size, mat_size, &a_pair)
        .expect("Failed to build sparse matrix pattern from triplets.");  // sparsity pattern
        let num_triplet = a_triplet.len();
        
        // convert triplets to matrix
        let a_vals: Vec<f64> = a_triplet.iter().map(|t| t.val).collect();
        let mut a_mat = SparseColMat::new_from_argsort(symbolic.clone(), &argsort, &a_vals)
            .expect("Failed to create sparse matrix from triplets.");  // initial matrix

        let time_1 = Instant::now();
        let time_initial = time_1.duration_since(time_0);

        // iterate to convergence
        let mut iter = 0;
        while iter < max_iter {
            let time_i0 = Instant::now();

            // solve A_k x_undamped = b_k for x_undamped
            // solve x_{k+1} = (1 - damp) * x_k + damp * x_undamped for x_damp
            x_udmp_vec = solver.solve(&a_mat, &b_vec, &x_iter_vec, mat_size);
            let x_damp_new = (1.0 - damp) * &x_iter_vec + damp * &x_udmp_vec;

            let time_i1 = Instant::now();

            // reassemble A_{k+1} and b_{k+1} with x_{k+1}
            a_triplet = Vec::with_capacity(num_triplet);  // reset A
            b_vec = Col::zeros(mat_size);  // reset b
            vars.update_unknown(&x_damp_new);
            self.assemble_matrix(vars, &mut a_triplet, &mut b_vec);
            let a_vals: Vec<f64> = a_triplet.iter().map(|t| t.val).collect();
            a_mat = SparseColMat::new_from_argsort(symbolic.clone(), &argsort, &a_vals)
                .expect("Failed to create sparse matrix from triplets.");

            let time_i2 = Instant::now();

            // compute residual
            let res = (&a_mat * &x_damp_new - &b_vec).norm_l2() / (b_vec.norm_l2() + 1e-10);
            println!("Iteration: {iter}; Residual: {res}");
            if res < tol {
                break;
            }
            x_iter_vec = x_damp_new;

            // update time measurements
            time_assemble += time_i2.duration_since(time_i1);
            time_solve += time_i1.duration_since(time_i0);

            // increment iteration
            iter += 1;
        }

        // error if not converged
        if iter == max_iter {
            return Err(FEChemError::FailedConvergence {
                caller: "SteadyBase::solve".to_string(),
                max_iter,
            });
        }

        let time_2 = Instant::now();

        // write variables
        vars.write_scalar(0.0, 0)?;

        let time_end = Instant::now();
        let time_write = time_end.duration_since(time_2);
        let time_total = time_end.duration_since(time_start);

        // output time measurement (Duration -> seconds as f64)
        println!("Solution completed!");
        println!("Total time: {:.6} s", time_total.as_secs_f64());
        println!("  Initialization time: {:.6} s", time_initial.as_secs_f64());
        println!("  Assembly time: {:.6} s", time_assemble.as_secs_f64());
        println!("  Solve time: {:.6} s", time_solve.as_secs_f64());
        println!("  Write time: {:.6} s", time_write.as_secs_f64());

        Ok(())
    }
}
