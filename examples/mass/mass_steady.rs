use fechem_fem2d::*;
use std::fs::create_dir_all;
use std::collections::HashMap;

/// Steady-state diffusion-reaction equation.
/// Run with: `cargo run --release --example mass_steady`
///
/// Geometry:
/// - x-axis: 0.0 m to 1.0 m; 20 elements
/// - y-axis: 0.0 m to 1.0 m; 20 elements
///
/// Properties:
/// - Diffusion coefficient (m2 s-1)
/// [D_AA D_AB] = [1.0 0.5]
/// [D_BA D_BB]   [0.0 2.0]
/// - Reaction rate (mol m-3 s-1; positive sign indicates production)
/// [R_A] = [-5.0 * c_A * c_B]
/// [R_B]   [+2.0 * c_A]
/// c_A and c_B in mol m-3
///
/// Boundary conditions:
/// - Left boundary (no flux): 0
/// - Bottom boundary (no flux): 0
/// - Right boundary (concentration):
/// [c_A] = [1.0 mol m-3]
/// [c_B]   [0.0 mol m-3]
/// - Top boundary (concentration):
/// [c_A] = [0.0 mol m-3]
/// [c_B]   [2.0 mol m-3]
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_mass_steady").unwrap();

    // problem and mesh
    // arguments: x_min, y_min, x_max, y_max, num_elem_x, num_elem_y
    let mut vars = Variables::new_from_bounds(0.0, 0.0, 1.0, 1.0, 20, 20)?;

    // geometry
    let dom = vars.add_dom(0)?;
    let bnd_l = vars.add_bnd(dom, 0)?;  // left
    let bnd_r = vars.add_bnd(dom, 1)?;  // right
    let bnd_b = vars.add_bnd(dom, 2)?;  // bottom
    let bnd_t = vars.add_bnd(dom, 3)?;  // top

    // variables
    // arguments: domain, initial_value, output_file
    // initial_value is an initial guess for steady-state problems
    let conc_a = vars.add_scldom_unk(dom, 0.0, "examples/output_mass_steady/conc_a.vtu".to_string())?;
    let conc_b = vars.add_scldom_unk(dom, 0.0, "examples/output_mass_steady/conc_b.vtu".to_string())?;

    // diffusion coefficients
    // no need to declare zero diffusion coefficients
    // must collect in hashmap with driving concentration index as key
    let diff_aa = vars.add_scldom_con(dom, 1.0, "".to_string())?;
    let diff_ab = vars.add_scldom_con(dom, 0.5, "".to_string())?;
    let diff_bb = vars.add_scldom_con(dom, 2.0, "".to_string())?;
    let diff_a = HashMap::from([(0, diff_aa), (1, diff_ab)]);
    let diff_b = HashMap::from([(1, diff_bb)]);

    // reaction rates (non-constant)
    // arguments: boundary, value_func, scldom_ids, output_file
    // value_func: time, scalar values (in the same order as scldom_ids; must be unknown-type scalars)
    let msrc_a_func = |_t: f64, scl:&[f64]| -5.0 * scl[0] * scl[1];
    let msrc_b_func = |_t: f64, scl:&[f64]| 2.0 * scl[0];
    let msrc_a = vars.add_scldom_fun(dom, Box::new(msrc_a_func), vec![conc_a, conc_b], "".to_string())?;
    let msrc_b = vars.add_scldom_fun(dom, Box::new(msrc_b_func), vec![conc_a], "".to_string())?;

    // constant properties on boundaries
    // arguments: boundary, value, output_file
    let n_a_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?; // no flux
    let n_b_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?; // no flux
    let c_a_r = vars.add_sclbnd_con(bnd_r, 1.0, "".to_string())?; // concentration
    let c_b_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?; // concentration
    let n_a_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?; // no flux
    let n_b_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?; // no flux
    let c_a_t = vars.add_sclbnd_con(bnd_t, 0.0, "".to_string())?; // concentration
    let c_b_t = vars.add_sclbnd_con(bnd_t, 2.0, "".to_string())?; // concentration

    // matrix solver
    // arguments: num_thread
    let solver = SolverLu::new(1)?;

    // physics solver
    let mut phys = SteadyMass::new(2);
    phys.add_mass_dom(0, dom, conc_a, diff_a, msrc_a);  // domain, component A
    phys.add_mass_dom(1, dom, conc_b, diff_b, msrc_b);  // domain, component B
    phys.add_conc_bnd(0, bnd_r, c_a_r);  // concentration BC, component A
    phys.add_conc_bnd(1, bnd_r, c_b_r);  // concentration BC, component B
    phys.add_conc_bnd(0, bnd_t, c_a_t);  // concentration BC, component A
    phys.add_conc_bnd(1, bnd_t, c_b_t);  // concentration BC, component B
    phys.add_mflx_bnd(0, bnd_l, n_a_l);  // flux BC, component A
    phys.add_mflx_bnd(1, bnd_l, n_b_l);  // flux BC, component B
    phys.add_mflx_bnd(0, bnd_b, n_a_b);  // flux BC, component A
    phys.add_mflx_bnd(1, bnd_b, n_b_b);  // flux BC, component B

    // solve
    // arguments: vars, solver, max_iter, tol, damping_factor
    // lower damping factor for better stability; higher for faster convergence
    phys.solve(&mut vars, Box::new(solver), 20, 1e-3, 0.8)?;

    Ok(())
}
