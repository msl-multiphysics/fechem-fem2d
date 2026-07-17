use fechem_fem2d::*;
use std::fs::create_dir_all;
use std::collections::HashMap;

/// Steady-state diffusion-reaction equation.
/// Run with: `cargo run --release --example mass_react`
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
/// [non-constant properties in SI units]
/// dom_0 - diffusion coefficients (m2 s-1)
/// [D_AA D_AB] = [1.0 0.5]
/// [D_BA D_BB]   [0.0 2.0]
/// dom_0 - reaction rates (mol m-3 s-1; positive if production)
/// [R_A] = [-5.0 * c_A * c_B]
/// [R_B]   [+2.0 * c_A]
///
/// Boundary conditions:
/// bnd_0 - no flux (mol m-2 s-1)
/// [N_A] = [0.0]
/// [N_B]   [0.0]
/// bnd_1 - concentration (mol m-3)
/// [c_A] = [1.0]
/// [c_B]   [0.0]
/// bnd_2 - no flux (mol m-2 s-1)
/// [N_A] = [0.0]
/// [N_B]   [0.0]
/// bnd_3 - concentration A; flux B (mol m-3; mol m-2 s-1)
/// [c_A] = [0.0]
/// [N_B]   [0.1]
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_mass_react").unwrap();

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
    let conc_a = vars.add_scldom_unk(dom, 0.0, "examples/output_mass_react/conc_a.vtu".to_string())?;
    let conc_b = vars.add_scldom_unk(dom, 0.0, "examples/output_mass_react/conc_b.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    // no need to declare zero diffusion coefficients
    // must collect in hashmap with driving concentration index as key
    let diff_aa = vars.add_scldom_con(dom, 1.0, "".to_string())?;  // D_AA
    let diff_ab = vars.add_scldom_con(dom, 0.5, "".to_string())?;  // D_AB
    let diff_bb = vars.add_scldom_con(dom, 2.0, "".to_string())?;  // D_BB
    let diff_a = HashMap::from([(0, diff_aa), (1, diff_ab)]);
    let diff_b = HashMap::from([(1, diff_bb)]);

    // non-constant domain scalars
    // arguments: domain, value_func, scalar_ids, output_file
    // value_func - returns value of scalar as a function of time and *unknown* scalars
    // - time is zero for steady-state problems
    // - scalar values are given in the same order as in scalar_ids
    let msrc_a_func = |_t: f64, scl: &[f64]| -5.0 * scl[0] * scl[1];  // scl[0] is c_A; scl[1] is c_B
    let msrc_b_func = |_t: f64, scl: &[f64]| 2.0 * scl[0];  // scl[0] is c_A
    let msrc_a = vars.add_scldom_fun(dom, Box::new(msrc_a_func), vec![conc_a, conc_b], "".to_string())?;  // reaction rate R_A
    let msrc_b = vars.add_scldom_fun(dom, Box::new(msrc_b_func), vec![conc_a], "".to_string())?;  // reaction rate R_B

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let mflx_a_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?;  // molar flux (positive if outward; negative if inward)
    let mflx_b_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?;  // molar flux
    let conc_a_r = vars.add_sclbnd_con(bnd_r, 1.0, "".to_string())?;  // concentration
    let conc_b_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?;  // concentration
    let mflx_a_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;  // molar flux
    let mflx_b_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;  // molar flux
    let conc_a_t = vars.add_sclbnd_con(bnd_t, 0.0, "".to_string())?;  // concentration
    let mflx_b_t = vars.add_sclbnd_con(bnd_t, 0.1, "".to_string())?;  // molar flux

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    let mut phys = SteadyMass::new(2);
    phys.add_mass_dom(0, dom, conc_a, diff_a, msrc_a);  // arguments: component, domain, conc, diff, msrc
    phys.add_mass_dom(1, dom, conc_b, diff_b, msrc_b);  // arguments: component, domain, conc, diff, msrc
    phys.add_mflx_bnd(0, bnd_l, mflx_a_l);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_l, mflx_b_l);  // arguments: component, boundary, mflx
    phys.add_conc_bnd(0, bnd_r, conc_a_r);  // arguments: component, boundary, conc
    phys.add_conc_bnd(1, bnd_r, conc_b_r);  // arguments: component, boundary, conc
    phys.add_mflx_bnd(0, bnd_b, mflx_a_b);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_b, mflx_b_b);  // arguments: component, boundary, mflx
    phys.add_conc_bnd(0, bnd_t, conc_a_t);  // arguments: component, boundary, conc
    phys.add_mflx_bnd(1, bnd_t, mflx_b_t);  // arguments: component, boundary, mflx

    // physics solver
    // arguments: max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 20, 1e-3, 0.8)?;

    Ok(())
}
