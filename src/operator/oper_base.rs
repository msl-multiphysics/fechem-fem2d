use crate::base::scl_dom::ScalarDomainType;
use crate::base::vars::Variables;
use faer::Col;
use faer::sparse::Triplet;

pub trait OperatorBase {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64);

    fn add_a(
        &self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        scl_row: usize, row: usize, scl_col: usize, col: usize, value: f64
    ) {
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

    fn add_b(
        &self, vars: &Variables, b_vec: &mut Col<f64>,
        scl_row: usize, row: usize, value: f64
    ) {
        let row_start = match vars.scl_dom[scl_row].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let xid_row = row + row_start;
        b_vec[xid_row] += value;
    }

}
