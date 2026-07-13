use fechem_fem2d::*;
use std::fs::create_dir_all;

/// Steady-state heat equation with multiple domains.
/// Run with: `cargo run --release --example heat_multi`
///
/// Geometry:
/// - Square with uniformly sized tri elements
/// - Has a larger diagonal going from the lower left to the upper right
/// - Has a smaller diagonal going from the lower right to the center
///
/// Properties:
/// - Thermal conductivity (all domains): 1.0 W m-1 K-1
/// - Heat source (all domains): -500.0 W m-3
///
/// Boundary conditions:
/// - Left boundary (outward flux): 50.0 W m-2
/// - Bottom boundary (no flux): 0
/// - Right boundary (temperature): 300 K
/// - Top boundary (temperature): 400 K
/// 
/// Interface conditions:
/// - Temperature and flux continuity at interfaces
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_heat_multi").unwrap();

    // problem and mesh
    let mut vars = Variables::new("examples/gmsh/gmsh_uniform_multi.msh".to_string())?;

    // geometry
    // for interfaces, the order of the domains does not matter
    let dom_tl = vars.add_dom(0)?;  // top-left triangle
    let dom_r = vars.add_dom(1)?;   // right triangle
    let dom_b = vars.add_dom(2)?;   // bottom triangle
    let bnd_l = vars.add_bnd(dom_tl, 0)?;  // left (attached to top-left triangle)
    let bnd_r = vars.add_bnd(dom_r, 1)?;   // right (attached to right triangle)
    let bnd_b = vars.add_bnd(dom_b, 2)?;   // bottom (attached to bottom triangle)
    let bnd_t = vars.add_bnd(dom_tl, 3)?;  // top (attached to top-left triangle)
    let itf_l1 = vars.add_itf(dom_tl, dom_b, 4)?;  // large diagonal (lower-left half; joins top-left and bottom)
    let itf_l2 = vars.add_itf(dom_tl, dom_r, 5)?;  // large diagonal (upper-right half; joins top-left and right)
    let itf_s = vars.add_itf(dom_b, dom_r, 6)?;  // small diagonal (joins bottom and right)

    // variables
    // lagrange multipliers are needed for continuity interfaces
    let temp_tl = vars.add_scldom_unk(dom_tl, 0.0, "examples/output_heat_multi/temp_tl.vtu".to_string())?;
    let temp_r = vars.add_scldom_unk(dom_r, 0.0, "examples/output_heat_multi/temp_r.vtu".to_string())?;
    let temp_b = vars.add_scldom_unk(dom_b, 0.0, "examples/output_heat_multi/temp_b.vtu".to_string())?;
    let lmd_l1 = vars.add_sclitf_unk(itf_l1, 0.0, "".to_string())?;  // lagrange multiplier for large diagonal (lower-left half)
    let lmd_l2 = vars.add_sclitf_unk(itf_l2, 0.0, "".to_string())?;  // lagrange multiplier for large diagonal (upper-right half)
    let lmd_s = vars.add_sclitf_unk(itf_s, 0.0, "".to_string())?;  // lagrange multiplier for small diagonal

    // constant properties on mesh
    // arguments: domain, value, output_file
    // set output file to empty string to not write to file
    let cond_tl = vars.add_scldom_con(dom_tl, 1.0, "".to_string())?;
    let cond_r = vars.add_scldom_con(dom_r, 1.0, "".to_string())?;
    let cond_b = vars.add_scldom_con(dom_b, 1.0, "".to_string())?;
    let hsrc_tl = vars.add_scldom_con(dom_tl, -500.0, "".to_string())?;
    let hsrc_r = vars.add_scldom_con(dom_r, -500.0, "".to_string())?;
    let hsrc_b = vars.add_scldom_con(dom_b, -500.0, "".to_string())?;

    // constant properties on boundaries
    // arguments: boundary, value, output_file
    let n_l = vars.add_sclbnd_con(bnd_l, 50.0, "".to_string())?; // positive for outward heat flux
    let t_r = vars.add_sclbnd_con(bnd_r, 300.0, "".to_string())?;
    let n_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?;
    let t_t = vars.add_sclbnd_con(bnd_t, 400.0, "".to_string())?;

    // matrix solver
    // arguments: num_thread
    let solver = SolverLu::new(1)?;

    // physics solver
    let mut phys = SteadyHeat::new();
    phys.add_heat_dom(dom_tl, temp_tl, cond_tl, hsrc_tl);
    phys.add_heat_dom(dom_r, temp_r, cond_r, hsrc_r);
    phys.add_heat_dom(dom_b, temp_b, cond_b, hsrc_b);
    phys.add_hflx_bnd(bnd_l, n_l);  // flux BC
    phys.add_temp_bnd(bnd_r, t_r);  // temperature BC
    phys.add_hflx_bnd(bnd_b, n_b);  // flux BC
    phys.add_temp_bnd(bnd_t, t_t);  // temperature BC
    phys.add_cont_itf(itf_l1, lmd_l1); // continuity
    phys.add_cont_itf(itf_l2, lmd_l2); // continuity
    phys.add_cont_itf(itf_s, lmd_s);   // continuity

    // solve
    // arguments: vars, solver, max_iter, tol, damping_factor
    // lower damping factor for better stability; higher for faster convergence
    phys.solve(&mut vars, Box::new(solver), 10, 1e-3, 1.0)?;

    Ok(())
}
