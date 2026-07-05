use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vars::Variables;
use faer::Col;
use faer::sparse::Triplet;

pub trait OperatorBase {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64);

    // scalar-scalar blocks
    fn add_a_sclscl(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, scl_row: usize, row: usize, scl_col: usize, col: usize, value: f64) {
        let row_start = match vars.scl_dom[scl_row].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let col_start = match vars.scl_dom[scl_col].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let xid_row = row + row_start;
        let xid_col = col + col_start;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }

    fn add_b_scl(&self, vars: &Variables, b_vec: &mut Col<f64>, scl_row: usize, row: usize, value: f64) {
        let row_start = match vars.scl_dom[scl_row].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let xid_row = row + row_start;
        b_vec[xid_row] += value;
    }

    // interface-interface blocks
    fn add_a_itfitf(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, itf_row: usize, row: usize, itf_col: usize, col: usize, value: f64) {
        let row_start = match vars.scl_itf[itf_row].scl_type {
            ScalarInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar interface type."),
        };
        let col_start = match vars.scl_itf[itf_col].scl_type {
            ScalarInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar interface type."),
        };
        let xid_row = row + row_start;
        let xid_col = col + col_start;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }

    fn add_b_itf(&self, vars: &Variables, b_vec: &mut Col<f64>, itf_row: usize, row: usize, value: f64) {
        let row_start = match vars.scl_itf[itf_row].scl_type {
            ScalarInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar interface type."),
        };
        let xid_row = row + row_start;
        b_vec[xid_row] += value;
    }

    // scalar-interface blocks
    fn add_a_itfscl(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, itf_row: usize, row: usize, scl_col: usize, col: usize, value: f64) {
        let row_start = match vars.scl_itf[itf_row].scl_type {
            ScalarInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar interface type."),
        };
        let col_start = match vars.scl_dom[scl_col].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let xid_row = row + row_start;
        let xid_col = col + col_start;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }

    fn add_a_sclitf(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, scl_row: usize, row: usize, itf_col: usize, col: usize, value: f64) {
        let row_start = match vars.scl_dom[scl_row].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let col_start = match vars.scl_itf[itf_col].scl_type {
            ScalarInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar interface type."),
        };
        let xid_row = row + row_start;
        let xid_col = col + col_start;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }
}
