use fechem_fem2d::*;
use std::fs::create_dir_all;
use std::collections::HashMap;

/// Steady-state heat, mass, and momentum transfer with lid-driven cavity flow.
/// Run with: `cargo run --release --example heatmassflow_lid`
///
/// Geometry:
/// Square domain with smaller tri elements near the walls.
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
/// dom_0 - density (rho = 1000.0 kg m-3)
/// dom_0 - heat capacity (cp = 0.1 J kg^-1 K^-1)
/// dom_0 - thermal conductivity (k = 0.1 W m-1 K-1)
/// dom_0 - heat source (Q = -100.0 W m-3)
/// 
/// Boundary conditions:
/// bnd_0 - temperature (T = 300 K)
/// bnd_1 - temperature (T = 300 K)
/// bnd_2 - temperature (T = 300 K)
/// bnd_3 - temperature (T = 400 K)
/// 
/// Mass Transfer
/// 
/// Properties:
/// [non-constant properties in SI units]
/// dom_0 - diffusion coefficients (m2 s-1)
/// [D_AA D_AB] = [1.0e-4 0.0]
/// [D_BA D_BB]   [0.0 1.0e-4]
/// dom_0 - Arrhenius reaction rates (mol m-3 s-1; positive if production)
/// k = A * exp(-Ea / (R * T))
/// A = 16.8 s-1; Ea = 3.0e4 J mol-1; R = 8.314 J mol-1 K-1
/// [R_A] = [-k * c_A]
/// [R_B]   [+k * c_A]
///
/// Boundary conditions:
/// bnd_0 - concentration (mol m-3)
/// [c_A] = [1.0]
/// [c_B]   [0.0]
/// bnd_1 - concentration (mol m-3)
/// [c_A] = [0.0]
/// [c_B]   [0.0]
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
/// dom_0 - density (same as in heat transfer)
/// dom_0 - viscosity (mu = 0.001 Pa s)
/// dom_0 - body force (f = <0.0, 0.0> N m-3)
///
/// Boundary conditions:
/// bnd_0 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_1 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_2 - no-slip (v = <0.0, 0.0> m s-1)
/// bnd_3 - velocity (v = <4.0e-4, 0.0> m s-1)
/// Zero reference pressure at arbitrary point
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heatmassflow_lid").unwrap();

    // mesh and variables
    // new - import mesh from gmsh file
    // arguments: input_file
    let mut vars = Variables::new("examples/gmsh/gmsh_lid.msh".to_string())?;

    // geometry
    // gmsh counts 1D and 2D physical groups from 1
    // subtract 1 to get FEChem domain and boundary indices
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
    let temp = vars.add_scldom_unk(dom, 300.0, "examples/output_heatmassflow_lid/temp.vtu".to_string())?;

    // constant domain scalars
    // arguments: domain, value, output_file
    let vlcp = vars.add_scldom_con(dom, 1000.0 * 0.1, "".to_string())?;  // volumetric heat capacity (rho * cp)
    let cond = vars.add_scldom_con(dom, 0.1, "".to_string())?;  // thermal conductivity
    let hsrc = vars.add_scldom_con(dom, -100.0, "".to_string())?;  // heat source (positive if source; negative if sink)

    // constant boundary scalars
    // arguments: boundary, value, output_file
    let temp_l = vars.add_sclbnd_con(bnd_l, 300.0, "".to_string())?;  // temperature
    let temp_r = vars.add_sclbnd_con(bnd_r, 300.0, "".to_string())?;  // temperature
    let temp_b = vars.add_sclbnd_con(bnd_b, 300.0, "".to_string())?;  // temperature
    let temp_t = vars.add_sclbnd_con(bnd_t, 400.0, "".to_string())?;  // temperature

    // mass transfer

    // unknown domain scalars
    let conc_a = vars.add_scldom_unk(dom, 0.5, "examples/output_heatmassflow_lid/conc_a.vtu".to_string())?;
    let conc_b = vars.add_scldom_unk(dom, 0.0, "examples/output_heatmassflow_lid/conc_b.vtu".to_string())?;

    // constant domain scalars
    // no need to declare zero diffusion coefficients
    // must collect in hashmap with driving concentration index as key
    let diff_aa = vars.add_scldom_con(dom, 1.0e-4, "".to_string())?;  // D_AA
    let diff_bb = vars.add_scldom_con(dom, 1.0e-4, "".to_string())?;  // D_BB
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
        -16.8 * (-3.0e4 / (8.314 * temp)).exp() * c_a
    };  // scl[0] is c_A; scl[1] is T
    let msrc_b_func = |_t: f64, scl: &[f64]| {
        let (c_a, temp) = (scl[0], scl[1]);
        16.8 * (-3.0e4 / (8.314 * temp)).exp() * c_a
    };  // scl[0] is c_A; scl[1] is T
    let msrc_a = vars.add_scldom_fun(dom, Box::new(msrc_a_func), vec![conc_a, temp], "".to_string())?;  // reaction rate R_A
    let msrc_b = vars.add_scldom_fun(dom, Box::new(msrc_b_func), vec![conc_a, temp], "".to_string())?;  // reaction rate R_B

    // constant boundary scalars
    let conc_a_l = vars.add_sclbnd_con(bnd_l, 1.0, "".to_string())?;  // concentration
    let conc_b_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?;  // concentration
    let conc_a_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?;  // concentration
    let conc_b_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?;  // concentration
    let mflx_a_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;  // molar flux (positive if outward; negative if inward)
    let mflx_b_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;  // molar flux
    let mflx_a_t = vars.add_sclbnd_con(bnd_t, 0.0, "".to_string())?;  // molar flux
    let mflx_b_t = vars.add_sclbnd_con(bnd_t, 0.0, "".to_string())?;  // molar flux

    // momentum transfer

    // unknown domain vectors
    // arguments: domain, initial_value_x, initial_value_y, output_file
    // initial_value - initial guess for steady-state problems; initial_value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    let vel = vars.add_vecdom_unk(dom, 0.0, 0.0, "examples/output_heatmassflow_lid/vel.vtu".to_string())?;

    // unknown domain scalars
    let pres = vars.add_scldom_unk(dom, 0.0, "examples/output_heatmassflow_lid/pres.vtu".to_string())?;

    // constant domain scalars
    let den = vars.add_scldom_con(dom, 1000.0, "".to_string())?;  // density
    let visc = vars.add_scldom_con(dom, 0.001, "".to_string())?;  // viscosity

    // constant domain vectors
    // arguments: domain, value_x, value_y, output_file
    let fce = vars.add_vecdom_con(dom, 0.0, 0.0, "".to_string())?;  // body force (den * g; not acceleration g)

    // constant boundary vectors
    // arguments: boundary, value_x, value_y, output_file
    let vel_l = vars.add_vecbnd_con(bnd_l, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_r = vars.add_vecbnd_con(bnd_r, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_b = vars.add_vecbnd_con(bnd_b, 0.0, 0.0, "".to_string())?;  // velocity
    let vel_t = vars.add_vecbnd_con(bnd_t, 4.0e-4, 0.0, "".to_string())?;  // velocity

    // steady-state heat-mass-momentum transfer solver
    // SteadyHeatMassFlow::new - arguments: number of components
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // set_pres_ref - set reference pressure at a point
    let mut phys = SteadyHeatMassFlow::new(2);
    phys.add_heat_dom(dom, temp, vlcp, cond, hsrc);  // arguments: domain, T, rho*cp, k, Q
    phys.add_temp_bnd(bnd_l, temp_l);  // arguments: boundary, T
    phys.add_temp_bnd(bnd_r, temp_r);  // arguments: boundary, T
    phys.add_temp_bnd(bnd_b, temp_b);  // arguments: boundary, T
    phys.add_temp_bnd(bnd_t, temp_t);  // arguments: boundary, T
    phys.add_mass_dom(0, dom, conc_a, diff_a, msrc_a);  // arguments: component, domain, conc, diff, msrc
    phys.add_mass_dom(1, dom, conc_b, diff_b, msrc_b);  // arguments: component, domain, conc, diff, msrc
    phys.add_conc_bnd(0, bnd_l, conc_a_l);  // arguments: component, boundary, conc
    phys.add_conc_bnd(1, bnd_l, conc_b_l);  // arguments: component, boundary, conc
    phys.add_conc_bnd(0, bnd_r, conc_a_r);  // arguments: component, boundary, conc
    phys.add_conc_bnd(1, bnd_r, conc_b_r);  // arguments: component, boundary, conc
    phys.add_mflx_bnd(0, bnd_b, mflx_a_b);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_b, mflx_b_b);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(0, bnd_t, mflx_a_t);  // arguments: component, boundary, mflx
    phys.add_mflx_bnd(1, bnd_t, mflx_b_t);  // arguments: component, boundary, mflx
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
    phys.solve(&mut vars, Box::new(linsolve), 40, 1e-3, 0.8)?;

    Ok(())
}
