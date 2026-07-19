use fechem_fem2d::*;
use std::fs::create_dir_all;
use std::collections::HashMap;

/// Steady-state heat, mass, and momentum transfer through a heated channel.
/// Run with: `cargo run --release --example heatmassflow_heater`
///
/// Geometry:
/// Channel with bump surrounded by solid regions.
///
/// Domains:
/// 0 - channel
/// 1 - bottom solid
/// 2 - top solid
/// 
/// Boundaries:
/// 0 - inlet
/// 1 - outlet
/// 2 - bottom of solid
/// 3 - top of solid
/// 
/// Interfaces:
/// 4 - bottom of channel
/// 5 - top of channel
///
/// Heat Transfer
/// 
/// Properties:
/// dom_0 - density (rho = 1000.0 kg m-3)
/// dom_0 - heat capacity (cp = 0.1 J kg^-1 K^-1)
/// dom_0 - thermal conductivity (k = 0.1 W m-1 K-1)
/// dom_0 - heat source (Q = 0.0 W m-3)
/// dom_1 - thermal conductivity (k = 1.0 W m-1 K-1)
/// dom_1 - heat source (Q = 50.0 W m-3)
/// dom_2 - thermal conductivity (k = 1.0 W m-1 K-1)
/// dom_2 - heat source (Q = 50.0 W m-3)
///
/// Boundary conditions:
/// bnd_0 - temperature (T = 300 K)
/// bnd_1 - outflow
/// bnd_2 - temperature (T = 400 K)
/// bnd_3 - temperature (T = 400 K)
///
/// Interface conditions:
/// itf_4 - contact resistance (0.1 W m-2 K-1)
/// itf_5 - contact resistance (0.1 W m-2 K-1)
/// 
/// Mass Transfer
/// 
/// Properties:
/// [non-constant properties in SI units]
/// dom_0 - diffusion coefficients (m2 s-1)
/// [D_AA D_AB] = [1.0e-3 0.0]
/// [D_BA D_BB]   [0.0 1.0e-3]
/// dom_0 - Arrhenius reaction rates (mol m-3 s-1; positive if production)
/// k = A * exp(-Ea / (R * T))
/// A = 168.0 s-1; Ea = 3.0e4 J mol-1; R = 8.314 J mol-1 K-1
/// [R_A] = [-k * c_A]
/// [R_B]   [+k * c_A]
/// dom_1 - unused domain; solid
/// dom_2 - unused domain; solid
///
/// Boundary conditions:
/// bnd_0 - concentration (mol m-3)
/// [c_A] = [1.0]
/// [c_B]   [0.0]
/// bnd_1 - outflow
/// bnd_4 - no flux (mol m-2 s-1)
/// [N_A] = [0.0]
/// [N_B]   [0.0]
/// bnd_5 - no flux (mol m-2 s-1)
/// [N_A] = [0.0]
/// [N_B]   [0.0]
/// 
/// Momentum Transfer
/// 
/// Properties:
/// dom_0 - density (same as in heat transfer)
/// dom_0 - viscosity (mu = 0.001 Pa s)
/// dom_0 - body force (f = <0.0, 0.0> N m-3)
/// dom_1 - unused domain; solid
/// dom_2 - unused domain; solid
///
/// Boundary conditions:
/// bnd_0 - velocity (v = <0.002, 0.0> m s-1)
/// bnd_1 - pressure (p = 0.0 Pa)
/// bnd_4 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_5 - no-slip (v = <0.0, 0.0> m s-1)
/// 
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heatmassflow_heater").unwrap();

    // mesh and variables
    // new - import mesh from gmsh file
    // arguments: input_file
    let mut vars = Variables::new("examples/gmsh/gmsh_heater.msh".to_string())?;

    // geometry
    // gmsh counts 1D and 2D physical groups from 1
    // subtract 1 to get FEChem domain and boundary indices
    // for interfaces, the order of the domains does not matter
    let dom_c = vars.add_dom(0)?;  // channel
    let dom_b = vars.add_dom(1)?;  // bottom solid
    let dom_t = vars.add_dom(2)?;  // top solid
    let bnd_i = vars.add_bnd(dom_c, 0)?;  // inlet
    let bnd_o = vars.add_bnd(dom_c, 1)?;  // outlet
    let bnd_bs = vars.add_bnd(dom_b, 2)?;  // bottom of solid
    let bnd_ts = vars.add_bnd(dom_t, 3)?;  // top of solid
    let itf_bc = vars.add_itf(dom_b, dom_c, 4)?;  // bottom of channel
    let itf_tc = vars.add_itf(dom_t, dom_c, 5)?;  // top of channel
    let bnd_bc = vars.add_bnd(dom_c, 4)?;  // channel bottom wall (no-slip)
    let bnd_tc = vars.add_bnd(dom_c, 5)?;  // channel top wall (no-slip)

    // heat transfer

    // unknown domain scalars
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let temp_c = vars.add_scldom_unk(dom_c, 300.0, "examples/output_heatmassflow_heater/temp_c.vtu".to_string())?;
    let temp_b = vars.add_scldom_unk(dom_b, 300.0, "examples/output_heatmassflow_heater/temp_b.vtu".to_string())?;
    let temp_t = vars.add_scldom_unk(dom_t, 300.0, "examples/output_heatmassflow_heater/temp_t.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    let vlcp_c = vars.add_scldom_con(dom_c, 1000.0 * 0.1, "".to_string())?;  // volumetric heat capacity (rho * cp)
    let vlcp_b = vars.add_scldom_con(dom_b, 1.0, "".to_string())?;  // unused for steady conduction-only solid
    let vlcp_t = vars.add_scldom_con(dom_t, 1.0, "".to_string())?;  // unused for steady conduction-only solid
    let cond_c = vars.add_scldom_con(dom_c, 0.1, "".to_string())?;  // thermal conductivity
    let cond_b = vars.add_scldom_con(dom_b, 1.0, "".to_string())?;  // thermal conductivity
    let cond_t = vars.add_scldom_con(dom_t, 1.0, "".to_string())?;  // thermal conductivity
    let hsrc_c = vars.add_scldom_con(dom_c, 0.0, "".to_string())?;  // heat source (positive if source; negative if sink)
    let hsrc_b = vars.add_scldom_con(dom_b, 50.0, "".to_string())?;  // heat source
    let hsrc_t = vars.add_scldom_con(dom_t, 50.0, "".to_string())?;  // heat source

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let temp_i = vars.add_sclbnd_con(bnd_i, 300.0, "".to_string())?;  // temperature
    let temp_bs = vars.add_sclbnd_con(bnd_bs, 400.0, "".to_string())?;  // temperature
    let temp_ts = vars.add_sclbnd_con(bnd_ts, 400.0, "".to_string())?;  // temperature

    // constant interface scalars
    // arguments: interface, value, output_file
    // contact resistance is needed for contact resistance interfaces
    let hres_bc = vars.add_sclitf_con(itf_bc, 0.1, "".to_string())?;  // contact resistance
    let hres_tc = vars.add_sclitf_con(itf_tc, 0.1, "".to_string())?;  // contact resistance

    // mass transfer

    // unknown domain scalars
    let conc_a = vars.add_scldom_unk(dom_c, 1.0, "examples/output_heatmassflow_heater/conc_a.vtu".to_string())?;
    let conc_b = vars.add_scldom_unk(dom_c, 0.0, "examples/output_heatmassflow_heater/conc_b.vtu".to_string())?;

    // constant domain scalars
    // no need to declare zero diffusion coefficients
    // must collect in hashmap with driving concentration index as key
    let diff_aa = vars.add_scldom_con(dom_c, 1.0e-3, "".to_string())?;  // D_AA
    let diff_bb = vars.add_scldom_con(dom_c, 1.0e-3, "".to_string())?;  // D_BB
    let diff_a = HashMap::from([(0, diff_aa)]);
    let diff_b = HashMap::from([(1, diff_bb)]);

    // non-constant domain scalars
    // arguments: domain, value_func, scalar_ids, output_file
    // value_func - returns value of scalar as a function of time and *unknown* scalars
    // - time is zero for steady-state problems
    // - scalar values are given in the same order as in scalar_ids
    // Arrhenius rate: k = A * exp(-Ea / (R * T)); R_A = -k * c_A; R_B = +k * c_A
    let msrc_a_func = |_t: f64, scl: &[f64]| {
        let (c_a, temp) = (scl[0], scl[1]);
        -168.0 * (-3.0e4 / (8.314 * temp)).exp() * c_a
    };  // scl[0] is c_A; scl[1] is T
    let msrc_b_func = |_t: f64, scl: &[f64]| {
        let (c_a, temp) = (scl[0], scl[1]);
        168.0 * (-3.0e4 / (8.314 * temp)).exp() * c_a
    };  // scl[0] is c_A; scl[1] is T
    let msrc_a = vars.add_scldom_fun(dom_c, Box::new(msrc_a_func), vec![conc_a, temp_c], "".to_string())?;  // reaction rate R_A
    let msrc_b = vars.add_scldom_fun(dom_c, Box::new(msrc_b_func), vec![conc_a, temp_c], "".to_string())?;  // reaction rate R_B

    // constant boundary scalars
    let conc_a_i = vars.add_sclbnd_con(bnd_i, 1.0, "".to_string())?;  // concentration
    let conc_b_i = vars.add_sclbnd_con(bnd_i, 0.0, "".to_string())?;  // concentration
    let mflx_a_bc = vars.add_sclbnd_con(bnd_bc, 0.0, "".to_string())?;  // molar flux (positive if outward; negative if inward)
    let mflx_b_bc = vars.add_sclbnd_con(bnd_bc, 0.0, "".to_string())?;  // molar flux
    let mflx_a_tc = vars.add_sclbnd_con(bnd_tc, 0.0, "".to_string())?;  // molar flux
    let mflx_b_tc = vars.add_sclbnd_con(bnd_tc, 0.0, "".to_string())?;  // molar flux

    // momentum transfer

    // unknown domain vectors
    // arguments: domain, initial_value_x, initial_value_y, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let vel = vars.add_vecdom_unk(dom_c, 0.0, 0.0, "examples/output_heatmassflow_heater/vel.vtu".to_string())?;

    // unknown domain scalars
    let pres = vars.add_scldom_unk(dom_c, 0.0, "examples/output_heatmassflow_heater/pres.vtu".to_string())?;

    // constant domain scalars
    let den = vars.add_scldom_con(dom_c, 1000.0, "".to_string())?;  // density
    let visc = vars.add_scldom_con(dom_c, 0.001, "".to_string())?;  // viscosity

    // constant domain vectors
    // arguments: domain, value_x, value_y, output_file
    let fce = vars.add_vecdom_con(dom_c, 0.0, 0.0, "".to_string())?;  // body force (den * g; not acceleration g)

    // constant boundary vectors
    // arguments: boundary, value_x, value_y, output_file
    let vel_i = vars.add_vecbnd_con(bnd_i, 0.002, 0.0, "".to_string())?;  // velocity
    let vel_bc = vars.add_vecbnd_con(bnd_bc, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_tc = vars.add_vecbnd_con(bnd_tc, 0.0, 0.0, "".to_string())?;  // velocity

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let pres_o = vars.add_sclbnd_con(bnd_o, 0.0, "".to_string())?;  // pressure

    // steady-state heat-mass-momentum transfer solver
    // SteadyHeatMassFlow::new - arguments: number of components
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hout_bnd - register boundary with heat outflow
    // add_hres_itf - register contact resistance interface
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // add_mout_bnd - register boundary with mass outflow
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = SteadyHeatMassFlow::new(2);
    phys.add_heat_dom(dom_c, temp_c, vlcp_c, cond_c, hsrc_c);  // arguments: domain, T, rho*cp, k, Q
    phys.add_heat_dom(dom_b, temp_b, vlcp_b, cond_b, hsrc_b);  // arguments: domain, T, rho*cp, k, Q
    phys.add_heat_dom(dom_t, temp_t, vlcp_t, cond_t, hsrc_t);  // arguments: domain, T, rho*cp, k, Q
    phys.add_temp_bnd(bnd_i, temp_i);  // arguments: boundary, T
    phys.add_hout_bnd(bnd_o);  // arguments: boundary
    phys.add_temp_bnd(bnd_bs, temp_bs);  // arguments: boundary, T
    phys.add_temp_bnd(bnd_ts, temp_ts);  // arguments: boundary, T
    phys.add_hres_itf(itf_bc, hres_bc);  // arguments: interface, contact resistance
    phys.add_hres_itf(itf_tc, hres_tc);  // arguments: interface, contact resistance
    phys.add_mass_dom(0, dom_c, conc_a, diff_a, msrc_a);  // arguments: component, domain, conc, diff, msrc
    phys.add_mass_dom(1, dom_c, conc_b, diff_b, msrc_b);  // arguments: component, domain, conc, diff, msrc
    phys.add_conc_bnd(0, bnd_i, conc_a_i);  // arguments: component, boundary, conc
    phys.add_conc_bnd(1, bnd_i, conc_b_i);  // arguments: component, boundary, conc
    phys.add_mout_bnd(0, bnd_o);  // arguments: component, boundary
    phys.add_mout_bnd(1, bnd_o);  // arguments: component, boundary
    phys.add_mflx_bnd(0, bnd_bc, mflx_a_bc);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_bc, mflx_b_bc);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(0, bnd_tc, mflx_a_tc);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_tc, mflx_b_tc);  // arguments: component, boundary, mflx
    phys.add_flow_dom(dom_c, vel, pres, den, visc, fce);  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i);  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o);  // arguments: boundary, p
    phys.add_vel_bnd(bnd_bc, vel_bc);  // arguments: boundary, v
    phys.add_vel_bnd(bnd_tc, vel_tc);  // arguments: boundary, v

    // physics solver
    // arguments: max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let linsolve = SolverLu::new(1)?;
    phys.solve(&mut vars, Box::new(linsolve), 20, 1e-5, 1.0)?;

    Ok(())
}
