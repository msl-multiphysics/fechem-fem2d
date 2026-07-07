use crate::base::scl_dom::ScalarDomainType;
use crate::base::scl_itf::ScalarInterfaceType;
use crate::base::vec_dom::VectorDomainType;
use crate::base::vec_itf::VectorInterfaceType;
use crate::base::vars::Variables;
use faer::Col;
use faer::sparse::Triplet;

pub trait OperatorBase {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, t: f64, factor: f64);

    // scalar domain-domain blocks
    fn add_a_scldom(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, scl_row: usize, row: usize, scl_col: usize, col: usize, value: f64) {
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
    fn add_b_scldom(&self, vars: &Variables, b_vec: &mut Col<f64>, scl_row: usize, row: usize, value: f64) {
        let row_start = match vars.scl_dom[scl_row].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let xid_row = row + row_start;
        b_vec[xid_row] += value;
    }

    // scalar interface-interface blocks
    fn add_a_sclitf(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, itf_row: usize, row: usize, itf_col: usize, col: usize, value: f64) {
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

    fn add_b_sclitf(&self, vars: &Variables, b_vec: &mut Col<f64>, itf_row: usize, row: usize, value: f64) {
        let row_start = match vars.scl_itf[itf_row].scl_type {
            ScalarInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar interface type."),
        };
        let xid_row = row + row_start;
        b_vec[xid_row] += value;
    }

    // vector domain-domain blocks
    fn add_a_vecdom(
        &self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        vec_row: usize, comp_row: usize, row: usize, vec_col: usize, comp_col: usize, col: usize, value: f64
    ) {
        let row_start = match vars.vec_dom[vec_row].vec_type {
            VectorDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown vector domain type."),
        };
        let col_start = match vars.vec_dom[vec_col].vec_type {
            VectorDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown vector domain type."),
        };
        let num_node = vars.dom[vec_row].num_node;
        let xid_row = row + row_start + comp_row * num_node;
        let xid_col = col + col_start + comp_col * num_node;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }
    fn add_b_vecdom(&self, vars: &Variables, b_vec: &mut Col<f64>, vec_row: usize, comp_row: usize, row: usize, value: f64) {
        let row_start = match vars.vec_dom[vec_row].vec_type {
            VectorDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown vector domain type."),
        };
        let num_node = vars.dom[vec_row].num_node;
        let xid_row = row + row_start + comp_row * num_node;
        b_vec[xid_row] += value;
    }

    // vector interface-interface blocks
    fn add_a_vecitf(
        &self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        itf_row: usize, comp_row: usize, row: usize, itf_col: usize, comp_col: usize, col: usize, value: f64
    ) {
        let row_start = match vars.vec_itf[itf_row].vec_type {
            VectorInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown vector interface type."),
        };
        let col_start = match vars.vec_itf[itf_col].vec_type {
            VectorInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown vector interface type."),
        };
        let num_node = vars.itf[itf_row].num_node;
        let xid_row = row + row_start + comp_row * num_node;
        let xid_col = col + col_start + comp_col * num_node;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }
    fn add_b_vecitf(&self, vars: &Variables, b_vec: &mut Col<f64>, itf_row: usize, comp_row: usize, row: usize, value: f64) {
        let row_start = match vars.vec_itf[itf_row].vec_type {
            VectorInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown vector interface type."),
        };
        let num_node = vars.itf[itf_row].num_node;
        let xid_row = row + row_start + comp_row * num_node;
        b_vec[xid_row] += value;
    }

    // scalar-vector blocks
    fn add_a_scldom_vecdom(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, scl_row: usize, row: usize, vec_col: usize, comp_col: usize, col: usize, value: f64) {
        let row_start = match vars.scl_dom[scl_row].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let col_start = match vars.vec_dom[vec_col].vec_type {
            VectorDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown vector domain type."),
        };
        let num_node = vars.dom[vec_col].num_node;
        let xid_row = row + row_start;
        let xid_col = col + col_start + comp_col * num_node;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }
    fn add_a_vecdom_scldom(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, vec_row: usize, comp_row: usize, row: usize, scl_col: usize, col: usize, value: f64) {
        let row_start = match vars.vec_dom[vec_row].vec_type {
            VectorDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown vector domain type."),
        };
        let col_start = match vars.scl_dom[scl_col].scl_type {
            ScalarDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown scalar domain type."),
        };
        let num_node = vars.dom[vec_row].num_node;
        let xid_row = row + row_start + comp_row * num_node;
        let xid_col = col + col_start;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }

    // scalar-interface blocks
    fn add_a_sclitf_scldom(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, itf_row: usize, row: usize, scl_col: usize, col: usize, value: f64) {
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
    fn add_a_scldom_sclitf(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, scl_row: usize, row: usize, itf_col: usize, col: usize, value: f64) {
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

    // vector-interface blocks
    fn add_a_vecitf_vecdom(
        &self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        itf_row: usize, comp_row: usize, row: usize, vec_col: usize, comp_col: usize, col: usize, value: f64
    ) {
        let row_start = match vars.vec_itf[itf_row].vec_type {
            VectorInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown vector interface type."),
        };
        let col_start = match vars.vec_dom[vec_col].vec_type {
            VectorDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown vector domain type."),
        };
        let num_node_row = vars.itf[itf_row].num_node;
        let num_node_col = vars.dom[vec_col].num_node;
        let xid_row = row + row_start + comp_row * num_node_row;
        let xid_col = col + col_start + comp_col * num_node_col;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }
    fn add_a_vecdom_vecitf(
        &self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        vec_row: usize, comp_row: usize, row: usize, itf_col: usize, comp_col: usize, col: usize, value: f64
    ) {
        let row_start = match vars.vec_dom[vec_row].vec_type {
            VectorDomainType::Unknown { start } => start,
            _ => panic!("Expected unknown vector domain type."),
        };
        let col_start = match vars.vec_itf[itf_col].vec_type {
            VectorInterfaceType::Unknown { start } => start,
            _ => panic!("Expected unknown vector interface type."),
        };
        let num_node_row = vars.dom[vec_row].num_node;
        let num_node_col = vars.itf[itf_col].num_node;
        let xid_row = row + row_start + comp_row * num_node_row;
        let xid_col = col + col_start + comp_col * num_node_col;
        a_triplet.push(Triplet::new(xid_row, xid_col, value));
    }

}
