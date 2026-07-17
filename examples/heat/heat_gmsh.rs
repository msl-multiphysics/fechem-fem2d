use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state heat equation with custom mesh.
/// Run with: `cargo run --release --example heat_gmsh`
///
/// Geometry:
/// Square domain with uniformly sized tri elements.
/// 
///          bnd_3
///       +---------+
///       |         |
/// bnd_0 |  dom_0  | bnd_1
///       |         |
///       +---------+
///          bnd_2
///
/// Domains:
/// 0 - main domain (1.0 m x 1.0 m)
/// 
/// Boundaries:
/// 0 - left boundary (1.0 m)
/// 1 - right boundary (1.0 m)
/// 2 - bottom boundary (1.0 m)
/// 3 - top boundary (1.0 m)
///
/// Properties:
/// dom_0 - thermal conductivity (k = 1.0 W m-1 K-1)
/// dom_0 - heat source (Q = -500.0 W m-3)
///
/// Boundary conditions:
/// bnd_0 - outward flux (q = 50.0 W m-2)
/// bnd_1 - temperature (T = 300 K)
/// bnd_2 - no flux (q = 0.0 W m-2)
/// bnd_3 - temperature (T = 400 K)
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heat_gmsh").unwrap();

    // mesh and variables
    // new - import mesh from gmsh file
    // arguments: input_file
    let mut vars = Variables::new("examples/gmsh/gmsh_uniform.msh".to_string())?;

    // geometry
    // gmsh counts 1D and 2D physical groups from 1
    // subtract 1 to get FEChem domain and boundary indices
    let dom = vars.add_dom(0)?;
    let bnd_l = vars.add_bnd(dom, 0)?;  // left
    let bnd_r = vars.add_bnd(dom, 1)?;  // right
    let bnd_b = vars.add_bnd(dom, 2)?;  // bottom
    let bnd_t = vars.add_bnd(dom, 3)?;  // top

    // unknown domain scalars
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let temp = vars.add_scldom_unk(dom, 0.0, "examples/output_heat_gmsh/temp.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    let cond = vars.add_scldom_con(dom, 1.0, "".to_string())?;  // thermal conductivity
    let hsrc = vars.add_scldom_con(dom, -500.0, "".to_string())?;  // heat source (positive if source; negative if sink)

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let hflx_l = vars.add_sclbnd_con(bnd_l, 50.0, "".to_string())?;  // heat flux (positive if outward; negative if inward)
    let temp_r = vars.add_sclbnd_con(bnd_r, 300.0, "".to_string())?;
    let hflx_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;
    let temp_t = vars.add_sclbnd_con(bnd_t, 400.0, "".to_string())?;

    // steady-state heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_hflx_bnd - register boundary with heat flux
    // add_temp_bnd - register boundary with temperature
    let mut phys = SteadyHeat::new();
    phys.add_heat_dom(dom, temp, cond, hsrc);  // arguments: domain, T, k, Q
    phys.add_hflx_bnd(bnd_l, hflx_l);  // arguments: boundary, q
    phys.add_temp_bnd(bnd_r, temp_r);  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_b, hflx_b);  // arguments: boundary, q
    phys.add_temp_bnd(bnd_t, temp_t);  // arguments: boundary, T

    // physics solver
    // arguments: max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 10, 1e-3, 1.0)?;

    Ok(())
}
