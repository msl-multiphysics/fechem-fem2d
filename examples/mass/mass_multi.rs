use fechem_fem2d::*;
use std::fs::create_dir_all;
use std::collections::HashMap;

/// Steady-state diffusion-reaction equation with multiple domains.
/// Run with: `cargo run --release --example mass_multi`
///
/// Geometry:
/// - Square with uniformly sized tri elements
/// - Has a larger diagonal going from the lower left to the upper right
/// - Has a smaller diagonal going from the lower right to the center
///
/// Properties:
/// - Diffusion coefficient (m2 s-1)
/// [D_AA D_AB] = [1.0 0.5]
/// [D_BA D_BB]   [0.0 2.0]
/// - Reaction rate (mol m-3 s-1; positive sign indicates production)
/// [R_A] = [-5.0 * c_A * c_B]
/// [R_B]   [+2.0 * c_A]
/// c_A and c_B in mol m-3
///
/// Boundary conditions:
/// - Left boundary (no flux): 0
/// - Bottom boundary (no flux): 0
/// - Right boundary (concentration):
/// [c_A] = [1.0 mol m-3]
/// [c_B]   [0.0 mol m-3]
/// - Top boundary (concentration):
/// [c_A] = [0.0 mol m-3]
/// [c_B]   [2.0 mol m-3]
///
/// Interface conditions:
/// - Concentration and flux continuity at interfaces
///
fn main() -> Result<(), FEChemError> {
    // output directory
    create_dir_all("examples/output_mass_multi").unwrap();

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
    // arguments: domain, initial_value, output_file
    // initial_value is an initial guess for steady-state problems
    // lagrange multipliers are needed for continuity interfaces
    let conc_a_tl = vars.add_scldom_unk(dom_tl, 0.0, "examples/output_mass_multi/conc_a_tl.vtu".to_string())?;
    let conc_b_tl = vars.add_scldom_unk(dom_tl, 0.0, "examples/output_mass_multi/conc_b_tl.vtu".to_string())?;
    let conc_a_r = vars.add_scldom_unk(dom_r, 0.0, "examples/output_mass_multi/conc_a_r.vtu".to_string())?;
    let conc_b_r = vars.add_scldom_unk(dom_r, 0.0, "examples/output_mass_multi/conc_b_r.vtu".to_string())?;
    let conc_a_b = vars.add_scldom_unk(dom_b, 0.0, "examples/output_mass_multi/conc_a_b.vtu".to_string())?;
    let conc_b_b = vars.add_scldom_unk(dom_b, 0.0, "examples/output_mass_multi/conc_b_b.vtu".to_string())?;
    let lmd_a_l1 = vars.add_sclitf_unk(itf_l1, 0.0, "".to_string())?;  // lagrange multiplier for large diagonal (lower-left half), component A
    let lmd_b_l1 = vars.add_sclitf_unk(itf_l1, 0.0, "".to_string())?;  // lagrange multiplier for large diagonal (lower-left half), component B
    let lmd_a_l2 = vars.add_sclitf_unk(itf_l2, 0.0, "".to_string())?;  // lagrange multiplier for large diagonal (upper-right half), component A
    let lmd_b_l2 = vars.add_sclitf_unk(itf_l2, 0.0, "".to_string())?;  // lagrange multiplier for large diagonal (upper-right half), component B
    let lmd_a_s = vars.add_sclitf_unk(itf_s, 0.0, "".to_string())?;  // lagrange multiplier for small diagonal, component A
    let lmd_b_s = vars.add_sclitf_unk(itf_s, 0.0, "".to_string())?;  // lagrange multiplier for small diagonal, component B

    // diffusion coefficients
    // no need to declare zero diffusion coefficients
    // must collect in hashmap with driving concentration index as key
    let diff_aa_tl = vars.add_scldom_con(dom_tl, 1.0, "".to_string())?;
    let diff_ab_tl = vars.add_scldom_con(dom_tl, 0.5, "".to_string())?;
    let diff_bb_tl = vars.add_scldom_con(dom_tl, 2.0, "".to_string())?;
    let diff_a_tl = HashMap::from([(0, diff_aa_tl), (1, diff_ab_tl)]);
    let diff_b_tl = HashMap::from([(1, diff_bb_tl)]);
    let diff_aa_r = vars.add_scldom_con(dom_r, 1.0, "".to_string())?;
    let diff_ab_r = vars.add_scldom_con(dom_r, 0.5, "".to_string())?;
    let diff_bb_r = vars.add_scldom_con(dom_r, 2.0, "".to_string())?;
    let diff_a_r = HashMap::from([(0, diff_aa_r), (1, diff_ab_r)]);
    let diff_b_r = HashMap::from([(1, diff_bb_r)]);
    let diff_aa_b = vars.add_scldom_con(dom_b, 1.0, "".to_string())?;
    let diff_ab_b = vars.add_scldom_con(dom_b, 0.5, "".to_string())?;
    let diff_bb_b = vars.add_scldom_con(dom_b, 2.0, "".to_string())?;
    let diff_a_b = HashMap::from([(0, diff_aa_b), (1, diff_ab_b)]);
    let diff_b_b = HashMap::from([(1, diff_bb_b)]);

    // reaction rates (non-constant)
    // arguments: domain, value_func, scldom_ids, output_file
    // value_func: time, scalar values (in the same order as scldom_ids; must be unknown-type scalars)
    let msrc_a_func = |_t: f64, scl:&[f64]| -5.0 * scl[0] * scl[1];
    let msrc_b_func = |_t: f64, scl:&[f64]| 2.0 * scl[0];
    let msrc_a_tl = vars.add_scldom_fun(dom_tl, Box::new(msrc_a_func), vec![conc_a_tl, conc_b_tl], "".to_string())?;
    let msrc_b_tl = vars.add_scldom_fun(dom_tl, Box::new(msrc_b_func), vec![conc_a_tl], "".to_string())?;
    let msrc_a_r = vars.add_scldom_fun(dom_r, Box::new(msrc_a_func), vec![conc_a_r, conc_b_r], "".to_string())?;
    let msrc_b_r = vars.add_scldom_fun(dom_r, Box::new(msrc_b_func), vec![conc_a_r], "".to_string())?;
    let msrc_a_b = vars.add_scldom_fun(dom_b, Box::new(msrc_a_func), vec![conc_a_b, conc_b_b], "".to_string())?;
    let msrc_b_b = vars.add_scldom_fun(dom_b, Box::new(msrc_b_func), vec![conc_a_b], "".to_string())?;

    // constant properties on boundaries
    // arguments: boundary, value, output_file
    let n_a_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?; // no flux
    let n_b_l = vars.add_sclbnd_con(bnd_l, 0.0, "".to_string())?; // no flux
    let c_a_r = vars.add_sclbnd_con(bnd_r, 1.0, "".to_string())?; // concentration
    let c_b_r = vars.add_sclbnd_con(bnd_r, 0.0, "".to_string())?; // concentration
    let n_a_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?; // no flux
    let n_b_b = vars.add_sclbnd_con(bnd_b, 0.0, "".to_string())?; // no flux
    let c_a_t = vars.add_sclbnd_con(bnd_t, 0.0, "".to_string())?; // concentration
    let c_b_t = vars.add_sclbnd_con(bnd_t, 2.0, "".to_string())?; // concentration

    // matrix solver
    // arguments: num_thread
    let solver = SolverLu::new(1)?;

    // physics solver
    let mut phys = SteadyMass::new(2);
    phys.add_mass_dom(0, dom_tl, conc_a_tl, diff_a_tl, msrc_a_tl);  // top-left, component A
    phys.add_mass_dom(1, dom_tl, conc_b_tl, diff_b_tl, msrc_b_tl);  // top-left, component B
    phys.add_mass_dom(0, dom_r, conc_a_r, diff_a_r, msrc_a_r);  // right, component A
    phys.add_mass_dom(1, dom_r, conc_b_r, diff_b_r, msrc_b_r);  // right, component B
    phys.add_mass_dom(0, dom_b, conc_a_b, diff_a_b, msrc_a_b);  // bottom, component A
    phys.add_mass_dom(1, dom_b, conc_b_b, diff_b_b, msrc_b_b);  // bottom, component B
    phys.add_conc_bnd(0, bnd_r, c_a_r);  // concentration BC, component A
    phys.add_conc_bnd(1, bnd_r, c_b_r);  // concentration BC, component B
    phys.add_conc_bnd(0, bnd_t, c_a_t);  // concentration BC, component A
    phys.add_conc_bnd(1, bnd_t, c_b_t);  // concentration BC, component B
    phys.add_mflx_bnd(0, bnd_l, n_a_l);  // flux BC, component A
    phys.add_mflx_bnd(1, bnd_l, n_b_l);  // flux BC, component B
    phys.add_mflx_bnd(0, bnd_b, n_a_b);  // flux BC, component A
    phys.add_mflx_bnd(1, bnd_b, n_b_b);  // flux BC, component B
    phys.add_cont_itf(0, itf_l1, lmd_a_l1); // continuity, component A
    phys.add_cont_itf(1, itf_l1, lmd_b_l1); // continuity, component B
    phys.add_cont_itf(0, itf_l2, lmd_a_l2); // continuity, component A
    phys.add_cont_itf(1, itf_l2, lmd_b_l2); // continuity, component B
    phys.add_cont_itf(0, itf_s, lmd_a_s);   // continuity, component A
    phys.add_cont_itf(1, itf_s, lmd_b_s);   // continuity, component B

    // solve
    // arguments: vars, solver, max_iter, tol, damping_factor
    // lower damping factor for better stability; higher for faster convergence
    phys.solve(&mut vars, Box::new(solver), 20, 1e-3, 0.8)?;

    Ok(())
}
