use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Transient Navier-Stokes equation with von Karman vortex street.
/// Run with: `cargo run --release --example flow_vortex`
///
/// Geometry:
/// - Channel with a circular cylinder for a von Karman vortex street
/// - Has smaller tri elements near the cylinder and in its downstream wake
/// - Channel is 2.2 m long and 0.41 m high
/// - Cylinder is 0.2 m from the left boundary and 0.2 m from the bottom boundary
/// - Cylinder radius is 0.05 m
/// - Cylinder is centered in the channel
///
/// Properties:
/// - Density: 1000.0 kg m-3
/// - Viscosity: 0.001 Pa s
///
/// Boundary conditions:
/// - Left boundary (velocity): <0.02, 0.0> m s-1
/// - Bottom boundary (velocity): <0.0, 0.0> m s-1
/// - Right boundary (pressure): 0.0 Pa
/// - Top boundary (velocity): <0.0, 0.0> m s-1
/// - Cylinder boundary (velocity): <0.0, 0.0> m s-1
/// 
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_flow_vortex").unwrap();

    // problem and mesh
    let mut vars = Variables::new("examples/gmsh/gmsh_vortex.msh".to_string())?;

    // geometry
    let dom = vars.add_dom(0)?;
    let bnd_l = vars.add_bnd(dom, 0)?;  // left
    let bnd_r = vars.add_bnd(dom, 1)?;  // right
    let bnd_b = vars.add_bnd(dom, 2)?;  // bottom
    let bnd_t = vars.add_bnd(dom, 3)?;  // top
    let bnd_c = vars.add_bnd(dom, 4)?;  // cylinder

    // variables
    // vector arguments: domain, initial_value_x, initial_value_y, output_file
    let vel = vars.add_vecdom_unk(dom, 0.0, 0.0, "examples/output_flow_vortex/vel.vtu".to_string())?;
    let pres = vars.add_scldom_unk(dom, 0.0, "examples/output_flow_vortex/pres.vtu".to_string())?;

    // constant properties on mesh
    // arguments: domain, value, output_file
    // set output file to empty string to not write to file
    let den = vars.add_scldom_con(dom, 1000.0, "".to_string())?;
    let visc = vars.add_scldom_con(dom, 0.1, "".to_string())?;
    let fce = vars.add_vecdom_con(dom, 0.0, 0.0, "".to_string())?;  // this is body force (den * g) not acceleration (g)

    // constant properties on boundaries
    // vector arguments: boundary, value_x, value_y, output_file
    let v_l = vars.add_vecbnd_con(bnd_l, 1.0, 0.0, "".to_string())?;
    let p_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?;
    let v_b = vars.add_vecbnd_con(bnd_b, 0.0, 0.0, "".to_string())?;
    let v_t = vars.add_vecbnd_con(bnd_t, 0.0, 0.0, "".to_string())?;
    let v_c = vars.add_vecbnd_con(bnd_c, 0.0, 0.0, "".to_string())?;

    // matrix solver
    // arguments: num_thread
    let solver = SolverLu::new(1)?;

    // physics solver
    let mut phys = TransientFlow::new();
    phys.add_flow_dom(dom, vel, pres, den, visc, fce);
    phys.add_vel_bnd(bnd_l, v_l);
    phys.add_pres_bnd(bnd_r, p_r);
    phys.add_vel_bnd(bnd_b, v_b);
    phys.add_vel_bnd(bnd_t, v_t);
    phys.add_vel_bnd(bnd_c, v_c);

    // solve
    // arguments: vars, solver, max_iter, tol, damping_factor
    // lower damping factor for better stability; higher for faster convergence
    phys.solve(&mut vars, Box::new(solver), 0.01, 100, 1, 100, 1e-3, 1.0)?;

    Ok(())
}
