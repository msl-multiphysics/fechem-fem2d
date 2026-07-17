use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Transient heat equation.
/// Run with: `cargo run --release --example heat_time`
///
/// Geometry:
/// Square domain defined by lengths along each axis.
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
/// - dom_0 - density (rho = 2.0 kg m-3)
/// - dom_0 - specific heat (cp = 5.0 J kg-1 K-1)
/// - dom_0 - thermal conductivity (k = 1.0 W m-1 K-1)
/// - dom_0 - heat source (Q = -500.0 W m-3)
///
/// Initial conditions:
/// - dom_0 - temperature (T = 300.0 K)
///
/// Boundary conditions:
/// - bnd_0 - outward flux (q = 50.0 W m-2)
/// - bnd_1 - temperature (T = 300.0 K)
/// - bnd_2 - no flux (q = 0.0 W m-2)
/// - bnd_3 - temperature (T = 400.0 K)
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heat_time").unwrap();

    // mesh and variables
    // new_from_bounds - create rectangular mesh from bounding box
    // arguments: x_min, y_min, x_max, y_max, num_elem_x, num_elem_y
    let mut vars = Variables::new_from_bounds(0.0, 0.0, 1.0, 1.0, 20, 20)?;

    // geometry
    // if using new_from_bounds, the domain is region 0
    // the left, right, bottom, and top boundaries are regions 0, 1, 2, and 3
    let dom = vars.add_dom(0)?;
    let bnd_l = vars.add_bnd(dom, 0)?;  // left
    let bnd_r = vars.add_bnd(dom, 1)?;  // right
    let bnd_b = vars.add_bnd(dom, 2)?;  // bottom
    let bnd_t = vars.add_bnd(dom, 3)?;  // top

    // unknown domain scalars
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // timestep is automatically appended to end of output file name
    let temp = vars.add_scldom_unk(dom, 300.0, "examples/output_heat_time/temp.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    let vlcp = vars.add_scldom_con(dom, 2.0 * 5.0, "".to_string())?;  // volumetric heat capacity (rho * cp)
    let cond = vars.add_scldom_con(dom, 1.0, "".to_string())?;
    let hsrc = vars.add_scldom_con(dom, -500.0, "".to_string())?;

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let hflx_l = vars.add_sclbnd_con(bnd_l, 50.0, "".to_string())?;
    let temp_r = vars.add_sclbnd_con(bnd_r, 300.0, "".to_string())?;
    let hflx_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;
    let temp_t = vars.add_sclbnd_con(bnd_t, 400.0, "".to_string())?;

    // transient heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_hflx_bnd - register boundary with heat flux
    // add_temp_bnd - register boundary with temperature
    let mut phys = TransientHeat::new();
    phys.add_heat_dom(dom, temp, vlcp, cond, hsrc);  // arguments: domain, T, rho * cp, k, Q
    phys.add_hflx_bnd(bnd_l, hflx_l);  // arguments: boundary, q
    phys.add_temp_bnd(bnd_r, temp_r);  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_b, hflx_b);  // arguments: boundary, q
    phys.add_temp_bnd(bnd_t, temp_t);  // arguments: boundary, T

    // physics solver
    // arguments: dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    // num_ts_write - writes an output file every num_ts_write time steps
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 0.1, 100, 1, 10, 1e-3,  1.0)?;

    Ok(())
}
