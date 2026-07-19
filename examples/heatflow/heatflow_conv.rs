use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state heat and momentum transfer with buoyancy-driven natural convection.
/// Run with: `cargo run --release --example heatflow_conv`
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
/// Heat Transfer
/// 
/// Properties:
/// [non-constant properties in SI units]
/// dom_0 - density (rho = 1000.0 * (1.0 - 0.001 * (T - 300.0)) kg m-3)
/// dom_0 - heat capacity (cp = 0.1 J kg^-1 K^-1)
/// dom_0 - thermal conductivity (k = 0.1 W m-1 K-1)
/// dom_0 - heat source (Q = 0.0 W m-3)
///
/// Boundary conditions:
/// bnd_0 - temperature (T = 310 K)
/// bnd_1 - temperature (T = 290 K)
/// bnd_2 - heat flux (q = 0.0 W m-2)
/// bnd_3 - heat flux (q = 0.0 W m-2)
///
/// Momentum Transfer
/// 
/// Properties:
/// [non-constant properties in SI units]
/// dom_0 - density (same as in heat transfer)
/// dom_0 - viscosity (mu = 10.0 Pa s)
/// dom_0 - body force (f = <0.0, -rho * 9.81>)
///
/// Boundary conditions:
/// bnd_0 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_1 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_2 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_3 - no-slip (v = <0.0, 0.0> m s-1)
/// Zero reference pressure at arbitrary point
/// 
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_conv").unwrap();

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

    // heat transfer

    // unknown domain scalars
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let temp = vars.add_scldom_unk(dom, 300.0, "examples/output_heatflow_conv/temp.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    let vlcp_func = |_t: f64, scl: &[f64]| 1000.0 * (1.0 - 0.001 * (scl[0] - 300.0)) * 0.1;  // scl[0] is T
    let vlcp = vars.add_scldom_fun(dom, Box::new(vlcp_func), vec![temp], "".to_string())?;  // volumetric heat capacity (rho * cp); scalar_ids[0] is also T
    let cond = vars.add_scldom_con(dom, 0.1, "".to_string())?;  // thermal conductivity
    let hsrc = vars.add_scldom_con(dom, 0.0, "".to_string())?;  // heat source (positive if source; negative if sink)

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let temp_l = vars.add_sclbnd_con(bnd_l, 310.0, "".to_string())?;  // temperature
    let temp_r = vars.add_sclbnd_con(bnd_r, 290.0, "".to_string())?;  // temperature
    let hflx_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;  // heat flux
    let hflx_t = vars.add_sclbnd_con(bnd_t, 0.0, "".to_string())?;  // heat flux

    // momentum transfer

    // unknown domain vectors
    // arguments: domain, initial_value_x, initial_value_y, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let vel = vars.add_vecdom_unk(dom, 0.0, 0.0, "examples/output_heatflow_conv/vel.vtu".to_string())?;

    // unknown domain scalars
    let pres = vars.add_scldom_unk(dom, 0.0, "examples/output_heatflow_conv/pres.vtu".to_string())?;

    // constant domain scalars
    let den_func = |_t: f64, scl: &[f64]| 1000.0 * (1.0 - 0.001 * (scl[0] - 300.0));
    let den = vars.add_scldom_fun(dom, Box::new(den_func), vec![temp], "".to_string())?;  // density; scalar_ids[0] is also T
    let visc = vars.add_scldom_con(dom, 10.0, "".to_string())?;  // viscosity

    // function domain vectors
    // arguments: domain, value_func, scalar_ids, output_file
    let fce_func = |_t: f64, scl: &[f64]| {
        let rho = 1000.0 * (1.0 - 0.001 * (scl[0] - 300.0));
        (0.0, -rho * 9.81)
    };
    let fce = vars.add_vecdom_fun(dom, Box::new(fce_func), vec![temp], "".to_string())?;  // body force (den * g; not acceleration g)

    // constant boundary vectors
    // arguments: boundary, value_x, value_y, output_file
    let vel_l = vars.add_vecbnd_con(bnd_l, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_r = vars.add_vecbnd_con(bnd_r, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_b = vars.add_vecbnd_con(bnd_b, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_t = vars.add_vecbnd_con(bnd_t, 0.0, 0.0, "".to_string())?;  // velocity

    // steady-state heat-momentum transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hflx_bnd - register boundary with heat flux
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // set_pres_ref - set reference pressure at a point
    let mut phys = SteadyHeatFlow::new();
    phys.add_heat_dom(dom, temp, vlcp, cond, hsrc);  // arguments: domain, T, rho*cp, k, Q
    phys.add_temp_bnd(bnd_l, temp_l);  // arguments: boundary, T
    phys.add_temp_bnd(bnd_r, temp_r);  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_b, hflx_b);  // arguments: boundary, q
    phys.add_hflx_bnd(bnd_t, hflx_t);  // arguments: boundary, q
    phys.add_flow_dom(dom, vel, pres, den, visc, fce);  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l);  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r);  // arguments: boundary, v
    phys.add_vel_bnd(bnd_b, vel_b);  // arguments: boundary, v
    phys.add_vel_bnd(bnd_t, vel_t);  // arguments: boundary, v
    phys.set_pres_ref(dom, 0, 0.0);  // arguments: domain, point_index, p

    // physics solver
    // arguments: max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 20, 1e-3, 0.8)?;

    Ok(())
}
