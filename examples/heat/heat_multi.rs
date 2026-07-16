use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state heat equation with multiple domains.
/// Run with: `cargo run --release --example heat_multi`
///
/// Geometry:
/// - Square with uniformly sized tri elements
/// - Has bottom-left and top-right corner rectangles
///
/// Properties:
/// - Thermal conductivity (middle): 0.5 W m-1 K-1
/// - Thermal conductivity (bottom): 1.0 W m-1 K-1
/// - Thermal conductivity (top): 1.0 W m-1 K-1
/// - Heat source (middle): 0.0 W m-3
/// - Heat source (bottom): -500.0 W m-3
/// - Heat source (top): +500.0 W m-3
///
/// Boundary conditions:
/// - External boundies (outward flux): 300 K
/// 
/// Interface conditions:
/// - Temperature and flux continuity at bottom-middle interface
/// - Contact resistance (0.1 W m-2 K-1) at top-middle interface
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heat_multi").unwrap();

    // problem and mesh
    let mut vars = Variables::new("examples/gmsh/gmsh_threereg.msh".to_string())?;

    // geometry
    // for interfaces, the order of the domains does not matter
    let dom_m = vars.add_dom(0)?;  // middle
    let dom_b = vars.add_dom(1)?;  // bottom
    let dom_t = vars.add_dom(2)?;  // top
    let bnd_lb = vars.add_bnd(dom_b, 0)?;  // bottom-left
    let bnd_br = vars.add_bnd(dom_m, 1)?;  // bottom-right
    let bnd_rt = vars.add_bnd(dom_t, 2)?;  // top-right
    let bnd_tb = vars.add_bnd(dom_m, 3)?;  // top-bottom
    let itf_bm = vars.add_itf(dom_b, dom_m, 4)?;  // bottom-middle
    let itf_tm = vars.add_itf(dom_t, dom_m, 5)?;  // top-middle

    // variables
    // lagrange multipliers are needed for continuity interfaces
    // however, they are not needed for contact resistance interfaces
    let temp_m = vars.add_scldom_unk(dom_m, 0.0, "examples/output_heat_multi/temp_m.vtu".to_string())?;
    let temp_b = vars.add_scldom_unk(dom_b, 0.0, "examples/output_heat_multi/temp_b.vtu".to_string())?;
    let temp_t = vars.add_scldom_unk(dom_t, 0.0, "examples/output_heat_multi/temp_t.vtu".to_string())?;
    let lmd_bm = vars.add_sclitf_unk(itf_bm, 0.0, "".to_string())?;  // lagrange multiplier for bottom-middle

    // constant properties on mesh
    // arguments: domain, value, output_file
    // set output file to empty string to not write to file
    let cond_m = vars.add_scldom_con(dom_m, 0.5, "".to_string())?;
    let cond_b = vars.add_scldom_con(dom_b, 1.0, "".to_string())?;
    let cond_t = vars.add_scldom_con(dom_t, 1.0, "".to_string())?;
    let hsrc_m = vars.add_scldom_con(dom_m, 0.0, "".to_string())?;
    let hsrc_b = vars.add_scldom_con(dom_b, -500.0, "".to_string())?;
    let hsrc_t = vars.add_scldom_con(dom_t, 500.0, "".to_string())?;

    // constant properties on boundaries and interfaces
    // arguments: boundary, value, output_file
    let t_lb = vars.add_sclbnd_con(bnd_lb, 300.0, "".to_string())?;
    let t_br = vars.add_sclbnd_con(bnd_br, 300.0, "".to_string())?;
    let t_rt = vars.add_sclbnd_con(bnd_rt, 300.0, "".to_string())?;
    let t_tb = vars.add_sclbnd_con(bnd_tb, 300.0, "".to_string())?;
    let r_tm = vars.add_sclitf_con(itf_tm, 0.1, "".to_string())?;  // contact resistance for top-middle

    // matrix solver
    // arguments: num_thread
    let solver = SolverLu::new(1)?;

    // physics solver
    let mut phys = SteadyHeat::new();
    phys.add_heat_dom(dom_m, temp_m, cond_m, hsrc_m);
    phys.add_heat_dom(dom_b, temp_b, cond_b, hsrc_b);
    phys.add_heat_dom(dom_t, temp_t, cond_t, hsrc_t);
    phys.add_temp_bnd(bnd_lb, t_lb);  // flux BC
    phys.add_temp_bnd(bnd_br, t_br);  // flux BC
    phys.add_temp_bnd(bnd_rt, t_rt);  // flux BC
    phys.add_temp_bnd(bnd_tb, t_tb);  // flux BC
    phys.add_cont_itf(itf_bm, lmd_bm); // continuity
    phys.add_hres_itf(itf_tm, r_tm); // contact resistance

    // solve
    // arguments: vars, solver, max_iter, tol, damping_factor
    // lower damping factor for better stability; higher for faster convergence
    phys.solve(&mut vars, Box::new(solver), 10, 1e-3, 1.0)?;

    Ok(())
}
