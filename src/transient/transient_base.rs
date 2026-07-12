use crate::base::error::FEChemError;
use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vars::Variables;
use crate::base::vec_dom::VectorDomainType;
use crate::base::vec_itf::VectorInterfaceType;
use crate::solver::solver_base::SolverBase;
use faer::Col;
use faer::sparse::SparseColMat;
use std::time::{Duration, Instant};

// base trait for transient solvers.
pub trait TransientBase {
    // to be implemented in specific solver
    fn assemble_operator(&mut self, vars: &mut Variables, mat_size: &mut usize); // must set mat_size
    fn assemble_matrix(&self, vars: &Variables, a_mat: &mut SparseColMat<usize, f64>, b_vec: &mut Col<f64>, mat_size: usize, t: f64, dt: f64); // must reset a_mat and b_vec

    fn solve(&mut self, vars: &mut Variables, solver: Box<dyn SolverBase>, dt: f64, num_ts: usize, num_ts_write: usize, max_iter: usize, tol: f64, damp: f64) -> Result<(), FEChemError> {
        let time_start = Instant::now();
        println!("Starting transient solver.");

        // error handling
        if dt <= 0.0 {
            return Err(FEChemError::InvalidTimeStep {caller: "TransientBase::solve".to_string(), dt});
        }
        if num_ts == 0 {
            return Err(FEChemError::InvalidNumTimeSteps {caller: "TransientBase::solve".to_string(), num_ts});
        }
        if num_ts_write == 0 {
            return Err(FEChemError::InvalidWriteFrequency {caller: "TransientBase::solve".to_string(), num_ts_write});
        }
        if max_iter == 0 {
            return Err(FEChemError::InvalidMaxIter {caller: "TransientBase::solve".to_string(), max_iter});
        }
        if tol <= 0.0 {
            return Err(FEChemError::InvalidTolerance {caller: "TransientBase::solve".to_string(), tol});
        }
        if damp <= 0.0 || damp > 1.0 {
            return Err(FEChemError::InvalidDamping {caller: "TransientBase::solve".to_string(), damp});
        }

        // initialize time measurement
        let mut time_assemble = Duration::ZERO;
        let mut time_solve = Duration::ZERO;
        let mut time_write = Duration::ZERO;

        // initialize operators
        // also compute the matrix size
        let mut mat_size = 0;
        self.assemble_operator(vars, &mut mat_size);

        // initialize solver vectors
        let mut a_mat: SparseColMat<usize, f64> = SparseColMat::try_new_from_triplets(0, 0, &[])
            .expect("Failed to create empty sparse matrix.");
        let mut b_vec: Col<f64> = Col::zeros(mat_size);
        let mut x_udmp_vec: Col<f64>;
        let mut x_iter_vec: Col<f64> = Col::zeros(mat_size);

        // load initial unknown values into iteration vector
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

        let time_w0 = Instant::now();

        // write initial condition
        vars.write_scalar(0.0, 0)?;

        let time_w1 = Instant::now();
        time_write += time_w1.duration_since(time_w0);

        // iterate over time steps
        for ts in 0..num_ts {
            let time_0 = Instant::now();

            // initial assembly of A and b
            let t_next = (ts + 1) as f64 * dt; // backward Euler; use next time step for function evaluation
            vars.update_function(t_next);
            vars.update_unknown(&x_iter_vec);
            self.assemble_matrix(vars, &mut a_mat, &mut b_vec, mat_size, t_next, dt);

            let time_1 = Instant::now();
            time_assemble += time_1.duration_since(time_0);

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
                // assumed that a_mat and b_vec are reset within assemble_matrix
                vars.update_unknown(&x_damp_new);
                vars.update_function(t_next);
                self.assemble_matrix(vars, &mut a_mat, &mut b_vec, mat_size, t_next, dt);

                let time_i2 = Instant::now();

                // compute residual
                let res = (&a_mat * &x_damp_new - &b_vec).norm_l2() / (b_vec.norm_l2() + 1e-10);
                println!("Timestep: {ts}; Iteration: {iter}; Residual: {res}");

                // commit current iterate before checking convergence so that
                // x_iter_vec stays in sync with `vars` (which already holds
                // x_damp_new via update_unknown above)
                x_iter_vec = x_damp_new;

                // update time measurements before the convergence check so
                // the converged iteration's solve and reassembly are counted
                time_assemble += time_i2.duration_since(time_i1);
                time_solve += time_i1.duration_since(time_i0);

                if res < tol {
                    break;
                }

                // increment iteration
                iter += 1;
            }

            // error if not converged
            if iter == max_iter {
                return Err(FEChemError::FailedConvergence {
                    caller: "TransientBase::solve".to_string(),
                    max_iter,
                });
            }

            let time_3 = Instant::now();

            // update for next time step
            vars.update_prev();

            let time_4 = Instant::now();
            time_assemble += time_4.duration_since(time_3);

            // write the converged state at t = (ts + 1) * dt as snapshot
            // index ts + 1. Indexing is offset by one from the loop variable
            // so that snapshot 0 is the initial condition and snapshot k > 0
            // is the state after k completed steps.
            if (ts + 1) % num_ts_write == 0 {
                vars.write_scalar(t_next, ts + 1)?;
                let time_5 = Instant::now();
                time_write += time_5.duration_since(time_4);
            }
        }

        let time_end = Instant::now();
        let time_total = time_end.duration_since(time_start);

        // output time measurement (Duration -> seconds as f64)
        println!("Solution completed!");
        println!("Total time: {:.6} s", time_total.as_secs_f64());
        println!("  Assembly time: {:.6} s", time_assemble.as_secs_f64());
        println!("  Solve time: {:.6} s", time_solve.as_secs_f64());
        println!("  Write time: {:.6} s", time_write.as_secs_f64());

        Ok(())
    }
}
