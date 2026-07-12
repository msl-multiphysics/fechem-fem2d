use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::scl_bnd::ScalarBoundary;
use crate::base::scl_dom::ScalarDomain;
use crate::base::vec_bnd::VectorBoundary;
use crate::base::vec_dom::VectorDomain;
use std::fs::File;
use std::io::Write;

const VTU_FLOAT_DECIMALS: usize = 6;

fn format_vtu_f64(value: f64) -> String {
    if value == 0.0 || value.abs() < 1e-10 {
        return "0.0".to_string();
    }
    format!("{:.prec$}", value, prec = VTU_FLOAT_DECIMALS)
}

pub fn write_scldom_vtu(dom: &Domain, scldom: &ScalarDomain, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.vtu", scldom.file_name, ts);
    let caller = "write_scldom_vtu";
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_scldom_vtu".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    let mut connectivity: Vec<usize> = Vec::new();
    let mut offsets: Vec<usize> = Vec::with_capacity(dom.num_elem);
    let mut cell_types: Vec<u8> = Vec::with_capacity(dom.num_elem);
    let mut offset = 0usize;
    for eid in 0..dom.num_elem {
        let nnode = dom.elem_node[eid];
        for loc in 0..nnode {
            connectivity.push(dom.elem_node_id[eid][loc]);
        }
        offset += nnode;
        offsets.push(offset);
        let vtk_cell_type = match nnode {
            3 => 5u8, // VTK_TRIANGLE
            4 => 9u8, // VTK_QUAD
            _ => 7u8, // VTK_POLYGON
        };
        cell_types.push(vtk_cell_type);
    }

    writeln!(file, "<?xml version=\"1.0\"?>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  <UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        dom.num_node, dom.num_elem
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for nid in 0..dom.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            format_vtu_f64(dom.node_x[nid]),
            format_vtu_f64(dom.node_y[nid])
        )
        .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for node_id in &connectivity {
        write!(file, "{} ", node_id).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for off in &offsets {
        write!(file, "{} ", off).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for cell_type in &cell_types {
        write!(file, "{} ", cell_type).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <PointData Scalars=\"value\">").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"value\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for value in &scldom.node_value {
        write!(file, "{} ", format_vtu_f64(*value)).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </PointData>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "    </Piece>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  </UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "</VTKFile>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    Ok(())
}

pub fn write_vecdom_vtu(dom: &Domain, vecdom: &VectorDomain, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.vtu", vecdom.file_name, ts);
    let caller = "write_vecdom_vtu";
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_vecdom_vtu".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    let mut connectivity: Vec<usize> = Vec::new();
    let mut offsets: Vec<usize> = Vec::with_capacity(dom.num_elem);
    let mut cell_types: Vec<u8> = Vec::with_capacity(dom.num_elem);
    let mut offset = 0usize;
    for eid in 0..dom.num_elem {
        let nnode = dom.elem_node[eid];
        for loc in 0..nnode {
            connectivity.push(dom.elem_node_id[eid][loc]);
        }
        offset += nnode;
        offsets.push(offset);
        let vtk_cell_type = match nnode {
            3 => 5u8, // VTK_TRIANGLE
            4 => 9u8, // VTK_QUAD
            _ => 7u8, // VTK_POLYGON
        };
        cell_types.push(vtk_cell_type);
    }

    writeln!(file, "<?xml version=\"1.0\"?>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  <UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        dom.num_node, dom.num_elem
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for nid in 0..dom.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            format_vtu_f64(dom.node_x[nid]),
            format_vtu_f64(dom.node_y[nid])
        )
        .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for node_id in &connectivity {
        write!(file, "{} ", node_id).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for off in &offsets {
        write!(file, "{} ", off).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for cell_type in &cell_types {
        write!(file, "{} ", cell_type).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <PointData Vectors=\"value\">").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"value\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for nid in 0..dom.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            format_vtu_f64(vecdom.node_value_x[nid]),
            format_vtu_f64(vecdom.node_value_y[nid])
        )
        .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </PointData>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "    </Piece>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  </UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "</VTKFile>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    Ok(())
}

pub fn write_sclbnd_vtu(bnd: &Boundary, sclbnd: &ScalarBoundary, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.vtu", sclbnd.file_name, ts);
    let caller = "write_sclbnd_vtu";
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_sclbnd_vtu".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    let mut connectivity: Vec<usize> = Vec::new();
    let mut offsets: Vec<usize> = Vec::with_capacity(bnd.num_elem);
    let mut cell_types: Vec<u8> = Vec::with_capacity(bnd.num_elem);
    let mut offset = 0usize;
    for eid in 0..bnd.num_elem {
        let nnode = bnd.elem_node[eid];
        for loc in 0..nnode {
            connectivity.push(bnd.elem_node_id[eid][loc]);
        }
        offset += nnode;
        offsets.push(offset);
        let vtk_cell_type = match nnode {
            2 => 3u8, // VTK_LINE
            _ => 4u8, // VTK_POLY_LINE
        };
        cell_types.push(vtk_cell_type);
    }

    writeln!(file, "<?xml version=\"1.0\"?>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  <UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        bnd.num_node, bnd.num_elem
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for nid in 0..bnd.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            format_vtu_f64(bnd.node_x[nid]),
            format_vtu_f64(bnd.node_y[nid])
        )
        .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for node_id in &connectivity {
        write!(file, "{} ", node_id).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for off in &offsets {
        write!(file, "{} ", off).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for cell_type in &cell_types {
        write!(file, "{} ", cell_type).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <PointData Scalars=\"value\">").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"value\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for value in &sclbnd.node_value {
        write!(file, "{} ", format_vtu_f64(*value)).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </PointData>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "    </Piece>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  </UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "</VTKFile>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    Ok(())
}

pub fn write_vecbnd_vtu(bnd: &Boundary, vecbnd: &VectorBoundary, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.vtu", vecbnd.file_name, ts);
    let caller = "write_vecbnd_vtu";
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(FEChemError::FileWriteError {
                caller: "write_vecbnd_vtu".to_string(),
                file_path: file_path.clone(),
            });
        }
    };

    let mut connectivity: Vec<usize> = Vec::new();
    let mut offsets: Vec<usize> = Vec::with_capacity(bnd.num_elem);
    let mut cell_types: Vec<u8> = Vec::with_capacity(bnd.num_elem);
    let mut offset = 0usize;
    for eid in 0..bnd.num_elem {
        let nnode = bnd.elem_node[eid];
        for loc in 0..nnode {
            connectivity.push(bnd.elem_node_id[eid][loc]);
        }
        offset += nnode;
        offsets.push(offset);
        let vtk_cell_type = match nnode {
            2 => 3u8, // VTK_LINE
            _ => 4u8, // VTK_POLY_LINE
        };
        cell_types.push(vtk_cell_type);
    }

    writeln!(file, "<?xml version=\"1.0\"?>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  <UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        bnd.num_node, bnd.num_elem
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for nid in 0..bnd.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            format_vtu_f64(bnd.node_x[nid]),
            format_vtu_f64(bnd.node_y[nid])
        )
        .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Points>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for node_id in &connectivity {
        write!(file, "{} ", node_id).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for off in &offsets {
        write!(file, "{} ", off).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(
        file,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    write!(file, "          ").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for cell_type in &cell_types {
        write!(file, "{} ", cell_type).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file).map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </Cells>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "      <PointData Vectors=\"value\">").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"value\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    for nid in 0..bnd.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            format_vtu_f64(vecbnd.node_value_x[nid]),
            format_vtu_f64(vecbnd.node_value_y[nid])
        )
        .map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    }
    writeln!(file, "        </DataArray>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "      </PointData>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    writeln!(file, "    </Piece>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "  </UnstructuredGrid>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;
    writeln!(file, "</VTKFile>").map_err(|_| FEChemError::FileWriteError { caller: caller.to_string(), file_path: file_path.clone(), })?;

    Ok(())
}
