use crate::base::error::FEChemError;
use crate::base::geom_bnd::Boundary;
use crate::base::geom_dom::Domain;
use crate::base::scl_bnd::ScalarBoundary;
use crate::base::scl_dom::ScalarDomain;
use std::fs::File;
use std::io::Write;

pub fn write_scldom_vtu(dom: &Domain, scldom: &ScalarDomain, ts: usize) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.vtu", scldom.file_name, ts);
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
        let nnode = dom.elem_node_num[eid];
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

    writeln!(file, "<?xml version=\"1.0\"?>").expect("Unable to write VTU file");
    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )
    .expect("Unable to write VTU file");
    writeln!(file, "  <UnstructuredGrid>").expect("Unable to write VTU file");
    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        dom.num_node, dom.num_elem
    )
    .expect("Unable to write VTU file");

    writeln!(file, "      <Points>").expect("Unable to write VTU file");
    writeln!(
        file,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    for nid in 0..dom.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            dom.node_x[nid], dom.node_y[nid]
        )
        .expect("Unable to write VTU file");
    }
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");
    writeln!(file, "      </Points>").expect("Unable to write VTU file");

    writeln!(file, "      <Cells>").expect("Unable to write VTU file");
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for node_id in &connectivity {
        write!(file, "{} ", node_id).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");

    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for off in &offsets {
        write!(file, "{} ", off).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");

    writeln!(
        file,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for cell_type in &cell_types {
        write!(file, "{} ", cell_type).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");
    writeln!(file, "      </Cells>").expect("Unable to write VTU file");

    writeln!(file, "      <PointData Scalars=\"value\">").expect("Unable to write VTU file");
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"value\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for value in &scldom.node_value {
        write!(file, "{} ", value).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");
    writeln!(file, "      </PointData>").expect("Unable to write VTU file");

    writeln!(file, "    </Piece>").expect("Unable to write VTU file");
    writeln!(file, "  </UnstructuredGrid>").expect("Unable to write VTU file");
    writeln!(file, "</VTKFile>").expect("Unable to write VTU file");

    Ok(())
}

pub fn write_sclbnd_vtu(
    bnd: &Boundary,
    sclbnd: &ScalarBoundary,
    ts: usize,
) -> Result<(), FEChemError> {
    let file_path = format!("{}_{}.vtu", sclbnd.file_name, ts);
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
        let nnode = bnd.elem_node_num[eid];
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

    writeln!(file, "<?xml version=\"1.0\"?>").expect("Unable to write VTU file");
    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )
    .expect("Unable to write VTU file");
    writeln!(file, "  <UnstructuredGrid>").expect("Unable to write VTU file");
    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        bnd.num_node, bnd.num_elem
    )
    .expect("Unable to write VTU file");

    writeln!(file, "      <Points>").expect("Unable to write VTU file");
    writeln!(
        file,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    for nid in 0..bnd.num_node {
        writeln!(
            file,
            "          {} {} 0.0",
            bnd.node_x[nid], bnd.node_y[nid]
        )
        .expect("Unable to write VTU file");
    }
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");
    writeln!(file, "      </Points>").expect("Unable to write VTU file");

    writeln!(file, "      <Cells>").expect("Unable to write VTU file");
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for node_id in &connectivity {
        write!(file, "{} ", node_id).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");

    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for off in &offsets {
        write!(file, "{} ", off).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");

    writeln!(
        file,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for cell_type in &cell_types {
        write!(file, "{} ", cell_type).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");
    writeln!(file, "      </Cells>").expect("Unable to write VTU file");

    writeln!(file, "      <PointData Scalars=\"value\">").expect("Unable to write VTU file");
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"value\" format=\"ascii\">"
    )
    .expect("Unable to write VTU file");
    write!(file, "          ").expect("Unable to write VTU file");
    for value in &sclbnd.node_value {
        write!(file, "{} ", value).expect("Unable to write VTU file");
    }
    writeln!(file).expect("Unable to write VTU file");
    writeln!(file, "        </DataArray>").expect("Unable to write VTU file");
    writeln!(file, "      </PointData>").expect("Unable to write VTU file");

    writeln!(file, "    </Piece>").expect("Unable to write VTU file");
    writeln!(file, "  </UnstructuredGrid>").expect("Unable to write VTU file");
    writeln!(file, "</VTKFile>").expect("Unable to write VTU file");

    Ok(())
}
