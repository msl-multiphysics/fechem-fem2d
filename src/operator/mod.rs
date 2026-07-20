// base trait
pub mod oper_base;
pub mod prelude;

// scalar boundary operators
pub mod opscl_bnd_neu;
pub mod opscl_bnd_dir;
pub mod opscl_bnd_trn;
pub mod opscl_bnd_div;
pub mod opscl_bnd_out;
pub mod opscl_bnd_out_unity;

// scalar domain operators
pub mod opscl_dom_adv;
pub mod opscl_dom_adv_unity;
pub mod opscl_dom_supg_steady;
pub mod opscl_dom_supg_steady_unity;
pub mod opscl_dom_supg_time;
pub mod opscl_dom_supg_time_unity;
pub mod opscl_dom_diff;
pub mod opscl_dom_src;
pub mod opscl_dom_time;
pub mod opscl_dom_time_unity;
pub mod opscl_dom_den_time;
pub mod opscl_dom_div;
pub mod opscl_dom_pspg_steady;
pub mod opscl_dom_pspg_time;

// scalar interface operators
pub mod opscl_itf_cont;
pub mod opscl_itf_trn;

// vector domain operators
pub mod opvec_dom_adv;
pub mod opvec_dom_pres;
pub mod opvec_dom_diff;
pub mod opvec_dom_src;
pub mod opvec_dom_supg_steady;
pub mod opvec_dom_supg_time;
pub mod opvec_dom_time;

// vector boundary operators
pub mod opvec_bnd_dir;
pub mod opvec_bnd_pres;

// vector interface operators
pub mod opvec_itf_cont;
