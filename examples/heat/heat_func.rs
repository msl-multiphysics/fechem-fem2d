use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state heat equation with non-constant properties.
/// Run with: `cargo run --release --example heat_func`
///
/// Geometry:
/// - x-axis: 0.0 m to 1.0 m; 20 elements
/// - y-axis: 0.0 m to 1.0 m; 20 elements
///
/// Properties:
/// - Thermal conductivity: (0.1 + 0.3 * T[K]) W m-1 K-1
/// - Heat source: -(200.0 + 0.5 * T[K]) W m-3
///
/// Boundary conditions:
/// - Left boundary (outward flux): (10.0 + 0.1 * T[K]) W m-2
/// - Bottom boundary (no flux): 0
/// - Right boundary (temperature): 300 K
/// - Top boundary (temperature): 400 K
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heat_func").unwrap();

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
    let temp = vars.add_scldom_unk(dom, 0.0, "examples/output_heat_func/temp.vtu".to_string())?;

    // non-constant properties on mesh
    // arguments: domain, value_func, scldom_ids, output_file
    // value_func: time, scalar values (in the same order as scldom_ids; must be unknown-type scalars)
    let cond_func = |_t: f64, scl:&[f64]| 0.1 + 0.3 * scl[0];
    let hsrc_func = |_t: f64, scl:&[f64]| -(200.0 + 0.5 * scl[0]);
    let cond = vars.add_scldom_fun(dom, Box::new(cond_func), vec![temp], "".to_string())?;
    let hsrc = vars.add_scldom_fun(dom, Box::new(hsrc_func), vec![temp], "".to_string())?;

    // non-constant properties on boundaries
    // arguments: boundary, value_func, scldom_ids, output_file
    // value_func: time, scalar values (in the same order as scldom_ids; must be unknown-type scalars)
    let flux_func = |_t: f64, scl:&[f64]| 10.0 + 0.1 * scl[0];
    let n_l = vars.add_sclbnd_fun(bnd_l, Box::new(flux_func), vec![temp], "".to_string())?; // positive for outward heat flux
    let t_r = vars.add_sclbnd_con(bnd_r, 300.0, "".to_string())?;
    let n_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;
    let t_t = vars.add_sclbnd_con(bnd_t, 400.0, "".to_string())?;

    // matrix solver
    // arguments: num_thread
    let solver = SolverLu::new(1)?;

    // physics solver
    let mut phys = SteadyHeat::new();
    phys.add_heat_dom(dom, temp, cond, hsrc);
    phys.add_hflx_bnd(bnd_l, n_l); // flux BC
    phys.add_temp_bnd(bnd_r, t_r); // temperature BC
    phys.add_hflx_bnd(bnd_b, n_b); // flux BC
    phys.add_temp_bnd(bnd_t, t_t); // temperature BC

    // solve
    // arguments: vars, solver, max_iter, tol, damping_factor
    // lower damping factor for better stability; higher for faster convergence
    phys.solve(&mut vars, Box::new(solver), 10, 1e-3, 1.0)?;

    Ok(())
}
