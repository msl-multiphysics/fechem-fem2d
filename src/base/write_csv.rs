use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::scl_bnd::ScalarBoundary;
use crate::base::scl_dom::ScalarDomain;
use crate::base::vec_bnd::VectorBoundary;
use crate::base::vec_dom::VectorDomain;
use std::fs::File;
use std::io::Write;

pub fn write_scldom_csv(dom: &Domain, scldom: &ScalarDomain, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.csv", scldom.file_name, ts);
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_scldom_csv".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    writeln!(file, "x,y,value").expect("Unable to write header");
    for nid in 0..dom.num_node {
        writeln!(
            file,
            "{},{},{}",
            dom.node_x[nid], dom.node_y[nid], scldom.node_value[nid]
        )
        .expect("Unable to write data");
    }

    Ok(())
}

pub fn write_vecdom_csv(dom: &Domain, vecdom: &VectorDomain, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.csv", vecdom.file_name, ts);
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_vecdom_csv".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    writeln!(file, "x,y,value_x,value_y").expect("Unable to write header");
    for nid in 0..dom.num_node {
        writeln!(
            file,
            "{},{},{},{}",
            dom.node_x[nid],
            dom.node_y[nid],
            vecdom.node_value_x[nid],
            vecdom.node_value_y[nid]
        )
        .expect("Unable to write data");
    }

    Ok(())
}

pub fn write_sclbnd_csv(bnd: &Boundary, sclbnd: &ScalarBoundary, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.csv", sclbnd.file_name, ts);
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_sclbnd_csv".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    writeln!(file, "x,y,value").expect("Unable to write header");
    for nid in 0..bnd.num_node {
        writeln!(
            file,
            "{},{},{}",
            bnd.node_x[nid], bnd.node_y[nid], sclbnd.node_value[nid]
        )
        .expect("Unable to write data");
    }

    Ok(())
}

pub fn write_vecbnd_csv(bnd: &Boundary, vecbnd: &VectorBoundary, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.csv", vecbnd.file_name, ts);
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_vecbnd_csv".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    writeln!(file, "x,y,value_x,value_y").expect("Unable to write header");
    for nid in 0..bnd.num_node {
        writeln!(
            file,
            "{},{},{},{}",
            bnd.node_x[nid],
            bnd.node_y[nid],
            vecbnd.node_value_x[nid],
            vecbnd.node_value_y[nid]
        )
        .expect("Unable to write data");
    }

    Ok(())
}
