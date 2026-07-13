use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state heat equation with custom mesh.
/// Run with: `cargo run --release --example heat_gmsh`
///
/// Geometry:
/// - Square with uniformly sized tri elements
///
/// Properties:
/// - Thermal conductivity: 1.0 W m-1 K-1
/// - Heat source: -500.0 W m-3
///
/// Boundary conditions:
/// - Left boundary (outward flux): 50.0 W m-2
/// - Bottom boundary (no flux): 0
/// - Right boundary (temperature): 300 K
/// - Top boundary (temperature): 400 K
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heat_gmsh").unwrap();

    // problem and mesh
    let mut vars = Variables::new("examples/gmsh/gmsh_uniform.msh".to_string())?;

    // geometry
    let dom = vars.add_dom(0)?;
    let bnd_l = vars.add_bnd(dom, 0)?;  // left
    let bnd_r = vars.add_bnd(dom, 1)?;  // right
    let bnd_b = vars.add_bnd(dom, 2)?;  // bottom
    let bnd_t = vars.add_bnd(dom, 3)?;  // top

    // variables
    // arguments: domain, initial_value, output_file
    // initial_value is an initial guess for steady-state problems
    let temp = vars.add_scldom_unk(dom, 0.0, "examples/output_heat_gmsh/temp.vtu".to_string())?;

    // constant properties on mesh
    // arguments: domain, value, output_file
    // set output file to empty string to not write to file
    let cond = vars.add_scldom_con(dom, 1.0, "".to_string())?;
    let hsrc = vars.add_scldom_con(dom, -500.0, "".to_string())?;

    // constant properties on boundaries
    // arguments: boundary, value, output_file
    let n_l = vars.add_sclbnd_con(bnd_l, 50.0, "".to_string())?; // positive for outward heat flux
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
