use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state Navier-Stokes equation through a Z-shaped channel.
/// Run with: `cargo run --release --example flow_channel`
///
/// Geometry:
/// Z-shaped channel rotated 45 degrees counterclockwise.
/// Tube width is 0.2 m; each centerline segment is 0.5 m.
///
/// Domains:
/// 0 - main domain
/// 
/// Boundaries:
/// 0 - inlet
/// 1 - outlet
/// 2 - left wall
/// 3 - right wall
///
/// Properties:
/// dom_0 - density (rho = 1000.0 kg m-3)
/// dom_0 - viscosity (mu = 1.0 Pa s)
/// dom_0 - body force (f = <0.0, 0.0> N m-3)
///
/// Boundary conditions:
/// bnd_0 - velocity (v = <0.1/sqrt(2), 0.1/sqrt(2)> m s-1)
/// bnd_1 - pressure (p = 0.0 Pa)
/// bnd_2 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_3 - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_flow_channel").unwrap();

    // mesh and variables
    // new - import mesh from gmsh file
    // arguments: input_file
    let mut vars = Variables::new("examples/gmsh/gmsh_channel.msh".to_string())?;

    // geometry
    // gmsh counts 1D and 2D physical groups from 1
    // subtract 1 to get FEChem domain and boundary indices
    let dom = vars.add_dom(0)?;
    let bnd_i = vars.add_bnd(dom, 0)?;  // inlet
    let bnd_o = vars.add_bnd(dom, 1)?;  // outlet
    let bnd_l = vars.add_bnd(dom, 2)?;  // left wall
    let bnd_r = vars.add_bnd(dom, 3)?;  // right wall

    // unknown domain vectors
    // arguments: domain, initial_value_x, initial_value_y, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let vel = vars.add_vecdom_unk(dom, 0.0, 0.0, "examples/output_flow_channel/vel.vtu".to_string())?;

    // unknown domain scalars
    // arguments: domain, initial_value, output_file
    let pres = vars.add_scldom_unk(dom, 0.0, "examples/output_flow_channel/pres.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    let den = vars.add_scldom_con(dom, 1000.0, "".to_string())?;  // density
    let visc = vars.add_scldom_con(dom, 1.0, "".to_string())?;  // viscosity

    // constant domain vectors
    // arguments: domain, value_x, value_y, output_file
    let fce = vars.add_vecdom_con(dom, 0.0, 0.0, "".to_string())?;  // body force (den * g; not acceleration g)

    // constant boundary vectors
    // arguments: boundary, value_x, value_y, output_file
    let vel_i = vars.add_vecbnd_con(bnd_i, 0.1/f64::sqrt(2.0), 0.1/f64::sqrt(2.0), "".to_string())?;  // velocity
    let vel_l = vars.add_vecbnd_con(bnd_l, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_r = vars.add_vecbnd_con(bnd_r, 0.0, 0.0, "".to_string())?;  // velocity

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let pres_o = vars.add_sclbnd_con(bnd_o, 0.0, "".to_string())?;  // pressure

    // steady-state flow solver
    // add_flow_dom - register domain with flow
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = SteadyFlow::new();
    phys.add_flow_dom(dom, vel, pres, den, visc, fce);  // arguments: domain, vel, pres, den, visc, fce
    phys.add_vel_bnd(bnd_i, vel_i);  // arguments: boundary, vel
    phys.add_pres_bnd(bnd_o, pres_o);  // arguments: boundary, pres
    phys.add_vel_bnd(bnd_l, vel_l);  // arguments: boundary, vel
    phys.add_vel_bnd(bnd_r, vel_r);  // arguments: boundary, vel

    // physics solver
    // arguments: max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 20, 1e-3, 1.0)?;

    Ok(())
}
