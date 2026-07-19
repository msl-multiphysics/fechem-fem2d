use fechem_fem2d::*;
use std::fs::create_dir_all;
use std::collections::HashMap;

/// Steady-state mass and momentum transfer through a Z-shaped channel.
/// Run with: `cargo run --release --example massflow_channel`
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
/// Mass Transfer
/// 
/// Properties:
/// [non-constant properties in SI units]
/// dom_0 - diffusion coefficients (m2 s-1)
/// [D_AA D_AB] = [1.0e-3 0.0]
/// [D_BA D_BB]   [0.0 1.0e-3]
/// dom_0 - reaction rates (mol m-3 s-1; positive if production)
/// [R_A] = [-0.01 * c_A]
/// [R_B]   [+0.01 * c_A]
///
/// Boundary conditions:
/// bnd_0 - concentration (mol m-3)
/// [c_A] = [1.0]
/// [c_B]   [0.0]
/// bnd_1 - outflow
/// bnd_2 - no flux (mol m-2 s-1)
/// [N_A] = [0.0]
/// [N_B]   [0.0]
/// bnd_3 - no flux (mol m-2 s-1)
/// [N_A] = [0.0]
/// [N_B]   [0.0]
/// 
/// Momentum Transfer
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
    create_dir_all("examples/output_massflow_channel").unwrap();

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

    // mass transfer

    // unknown domain scalars
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let conc_a = vars.add_scldom_unk(dom, 1.0, "examples/output_massflow_channel/conc_a.vtu".to_string())?;
    let conc_b = vars.add_scldom_unk(dom, 0.0, "examples/output_massflow_channel/conc_b.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    // no need to declare zero diffusion coefficients
    // must collect in hashmap with driving concentration index as key
    let diff_aa = vars.add_scldom_con(dom, 1.0e-3, "".to_string())?;  // D_AA
    let diff_bb = vars.add_scldom_con(dom, 1.0e-3, "".to_string())?;  // D_BB
    let diff_a = HashMap::from([(0, diff_aa)]);
    let diff_b = HashMap::from([(1, diff_bb)]);

    // non-constant domain scalars
    // arguments: domain, value_func, scalar_ids, output_file
    // value_func - returns value of scalar as a function of time and *unknown* scalars
    // - time is zero for steady-state problems
    // - scalar values are given in the same order as in scalar_ids
    let msrc_a_func = |_t: f64, scl: &[f64]| -0.01 * scl[0];  // scl[0] is c_A
    let msrc_b_func = |_t: f64, scl: &[f64]| 0.01 * scl[0];  // scl[0] is c_A
    let msrc_a = vars.add_scldom_fun(dom, Box::new(msrc_a_func), vec![conc_a], "".to_string())?;  // reaction rate R_A
    let msrc_b = vars.add_scldom_fun(dom, Box::new(msrc_b_func), vec![conc_a], "".to_string())?;  // reaction rate R_B

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let conc_a_i = vars.add_sclbnd_con(bnd_i, 1.0, "".to_string())?;  // concentration
    let conc_b_i = vars.add_sclbnd_con(bnd_i, 0.0, "".to_string())?;  // concentration
    let mflx_a_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?;  // molar flux (positive if outward; negative if inward)
    let mflx_b_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?;  // molar flux
    let mflx_a_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?;  // molar flux
    let mflx_b_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?;  // molar flux

    // momentum transfer

    // unknown domain vectors
    // arguments: domain, initial_value_x, initial_value_y, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let vel = vars.add_vecdom_unk(dom, 0.0, 0.0, "examples/output_massflow_channel/vel.vtu".to_string())?;

    // unknown domain scalars
    let pres = vars.add_scldom_unk(dom, 0.0, "examples/output_massflow_channel/pres.vtu".to_string())?;

    // constant domain scalars
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

    // steady-state mass-momentum transfer solver
    // SteadyMassFlow::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // add_mout_bnd - register boundary with mass outflow
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = SteadyMassFlow::new(2);
    phys.add_mass_dom(0, dom, conc_a, diff_a, msrc_a);  // arguments: component, domain, conc, diff, msrc
    phys.add_mass_dom(1, dom, conc_b, diff_b, msrc_b);  // arguments: component, domain, conc, diff, msrc
    phys.add_conc_bnd(0, bnd_i, conc_a_i);  // arguments: component, boundary, conc
    phys.add_conc_bnd(1, bnd_i, conc_b_i);  // arguments: component, boundary, conc
    phys.add_mout_bnd(0, bnd_o);  // arguments: component, boundary
    phys.add_mout_bnd(1, bnd_o);  // arguments: component, boundary
    phys.add_mflx_bnd(0, bnd_l, mflx_a_l);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_l, mflx_b_l);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(0, bnd_r, mflx_a_r);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_r, mflx_b_r);  // arguments: component, boundary, mflx
    phys.add_flow_dom(dom, vel, pres, den, visc, fce);  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i);  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o);  // arguments: boundary, p
    phys.add_vel_bnd(bnd_l, vel_l);  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r);  // arguments: boundary, v

    // physics solver
    // arguments: max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 40, 1e-3, 0.8)?;

    Ok(())
}
