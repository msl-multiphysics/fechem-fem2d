use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state Navier-Stokes equation through a Z-shaped channel.
/// Run with: `cargo run --release --example flow_channel`
///
/// Geometry:
/// - Z-shaped channel
///
/// Properties:
/// - Density: 1000.0 kg m-3
/// - Viscosity: 1.0 Pa s
///
/// Boundary conditions:
/// - Inlet (velocity): <0.1/sqrt(2), 0.1/sqrt(2)> m s-1
/// - Right boundary (pressure): 0.0 Pa
/// 
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_flow_channel").unwrap();

    // problem and mesh
    let mut vars = Variables::new("examples/gmsh/gmsh_channel.msh".to_string())?;

    // geometry
    let dom = vars.add_dom(0)?;
    let bnd_i = vars.add_bnd(dom, 0)?;  // inlet
    let bnd_o = vars.add_bnd(dom, 1)?;  // outlet
    let bnd_l = vars.add_bnd(dom, 2)?;  // left
    let bnd_r = vars.add_bnd(dom, 3)?;  // right

    // variables
    // vector arguments: domain, initial_value_x, initial_value_y, output_file
    let vel = vars.add_vecdom_unk(dom, 0.0, 0.0, "examples/output_flow_channel/vel.vtu".to_string())?;
    let pres = vars.add_scldom_unk(dom, 0.0, "examples/output_flow_channel/pres.vtu".to_string())?;

    // constant properties on mesh
    // arguments: domain, value, output_file
    // set output file to empty string to not write to file
    let den = vars.add_scldom_con(dom, 1000.0, "".to_string())?;
    let visc = vars.add_scldom_con(dom, 1.0, "".to_string())?;
    let fce = vars.add_vecdom_con(dom, 0.0, 0.0, "".to_string())?;  // this is body force (den * g) not acceleration (g)

    // constant properties on boundaries
    // vector arguments: boundary, value_x, value_y, output_file
    let v_i = vars.add_vecbnd_con(bnd_i, 0.1/f64::sqrt(2.0), 0.1/f64::sqrt(2.0), "".to_string())?;
    let p_o = vars.add_sclbnd_con(bnd_o, 0.0, "".to_string())?;
    let v_l = vars.add_vecbnd_con(bnd_l, 0.0, 0.0, "".to_string())?;
    let v_r = vars.add_vecbnd_con(bnd_r, 0.0, 0.0, "".to_string())?;

    // matrix solver
    // arguments: num_thread
    let solver = SolverLu::new(1)?;

    // physics solver
    let mut phys = SteadyFlow::new();
    phys.add_flow_dom(dom, vel, pres, den, visc, fce);
    phys.add_vel_bnd(bnd_i, v_i);
    phys.add_pres_bnd(bnd_o, p_o);
    phys.add_vel_bnd(bnd_l, v_l);
    phys.add_vel_bnd(bnd_r, v_r);

    // solve
    // arguments: vars, solver, dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    // lower damping factor for better stability; higher for faster convergence
    // writes an output file every num_ts_write time steps
    phys.solve(&mut vars, Box::new(solver), 20, 1e-3, 1.0)?;

    Ok(())
}
