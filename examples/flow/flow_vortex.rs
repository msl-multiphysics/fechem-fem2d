use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Transient Navier-Stokes equation with von Karman vortex street.
/// Run with: `cargo run --release --example flow_vortex`
///
/// Geometry:
/// Channel with a circular cylinder; smaller tri elements near the
/// cylinder and in its downstream wake.
/// 
///                    bnd_3
///       +-----------------------------+
///       |                             |
/// bnd_0 |      (bnd_4)                | bnd_1
///       |                             |
///       +-----------------------------+
///                    bnd_2
///
/// Domains:
/// 0 - main domain (2.2 m x 0.41 m)
/// 
/// Boundaries:
/// 0 - left boundary (0.41 m)
/// 1 - right boundary (0.41 m)
/// 2 - bottom boundary (2.2 m)
/// 3 - top boundary (2.2 m)
/// 4 - cylinder boundary (radius = 0.05 m; center at <0.2, 0.2> m)
///
/// Properties:
/// dom_0 - density (rho = 1000.0 kg m-3)
/// dom_0 - viscosity (mu = 0.1 Pa s)
/// dom_0 - body force (f = <0.0, 0.0> N m-3)
///
/// Initial conditions:
/// dom_0 - velocity (v = <0.0, 0.0> m s-1)
/// dom_0 - pressure (p = 0.0 Pa)
///
/// Boundary conditions:
/// bnd_0 - velocity (v = <1.0, 0.0> m s-1)
/// bnd_1 - pressure (p = 0.0 Pa)
/// bnd_2 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_3 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_4 - no-slip (v = <0.0, 0.0> m s-1)
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_flow_vortex").unwrap();

    // mesh and variables
    // new - import mesh from gmsh file
    // arguments: input_file
    let mut vars = Variables::new("examples/gmsh/gmsh_vortex.msh".to_string())?;

    // geometry
    // gmsh counts 1D and 2D physical groups from 1
    // subtract 1 to get FEChem domain and boundary indices
    let dom = vars.add_dom(0)?;
    let bnd_l = vars.add_bnd(dom, 0)?;  // left
    let bnd_r = vars.add_bnd(dom, 1)?;  // right
    let bnd_b = vars.add_bnd(dom, 2)?;  // bottom
    let bnd_t = vars.add_bnd(dom, 3)?;  // top
    let bnd_c = vars.add_bnd(dom, 4)?;  // cylinder

    // unknown domain vectors
    // arguments: domain, initial_value_x, initial_value_y, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // timestep is automatically appended to end of output file name
    let vel = vars.add_vecdom_unk(dom, 0.0, 0.0, "examples/output_flow_vortex/vel.vtu".to_string())?;

    // unknown domain scalars
    // arguments: domain, initial_value, output_file
    let pres = vars.add_scldom_unk(dom, 0.0, "examples/output_flow_vortex/pres.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    let den = vars.add_scldom_con(dom, 1000.0, "".to_string())?;  // density
    let visc = vars.add_scldom_con(dom, 0.1, "".to_string())?;  // viscosity

    // constant domain vectors
    // arguments: domain, value_x, value_y, output_file
    let fce = vars.add_vecdom_con(dom, 0.0, 0.0, "".to_string())?;  // body force (den * g; not acceleration g)

    // constant boundary vectors
    // arguments: boundary, value_x, value_y, output_file
    let vel_l = vars.add_vecbnd_con(bnd_l, 1.0, 0.0, "".to_string())?;  // velocity
    let vel_b = vars.add_vecbnd_con(bnd_b, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_t = vars.add_vecbnd_con(bnd_t, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_c = vars.add_vecbnd_con(bnd_c, 0.0, 0.0, "".to_string())?;  // velocity

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let pres_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?;  // pressure

    // transient momentum transfer solver
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = TransientFlow::new();
    phys.add_flow_dom(dom, vel, pres, den, visc, fce);  // arguments: domain, vel, pres, den, visc, fce
    phys.add_vel_bnd(bnd_l, vel_l);  // arguments: boundary, vel
    phys.add_pres_bnd(bnd_r, pres_r);  // arguments: boundary, pres
    phys.add_vel_bnd(bnd_b, vel_b);  // arguments: boundary, vel
    phys.add_vel_bnd(bnd_t, vel_t);  // arguments: boundary, vel
    phys.add_vel_bnd(bnd_c, vel_c);  // arguments: boundary, vel

    // physics solver
    // arguments: dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    // num_ts_write - writes an output file every num_ts_write time steps
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 0.01, 100, 1, 100, 1e-3, 1.0)?;

    Ok(())
}
