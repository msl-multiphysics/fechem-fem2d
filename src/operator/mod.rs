// base trait
pub mod oper_base;
pub mod prelude;

// scalar boundary operators
pub mod opscl_bnd_neu;
pub mod opscl_bnd_dir;

// scalar domain operators
pub mod opscl_dom_diff;
pub mod opscl_dom_src;
pub mod opscl_dom_time;
// pub mod opscl_dom_div;
// pub mod opscl_dom_pspg;

// scalar interface operators
pub mod opscl_itf_cont;
// 
// // vector domain operators
// pub mod opvec_dom_adv;
// pub mod opvec_dom_pres;
// pub mod opvec_dom_diff;
// pub mod opvec_dom_src;
// pub mod opvec_dom_supg;
// 
// // vector boundary operators
// pub mod opvec_bnd_dir;
// 