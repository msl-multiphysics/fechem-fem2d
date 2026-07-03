//! Gmsh ASCII mesh reader for **format 4.1** (linear 2D elements only).
//! Supported element types: 1 (line), 2 (triangle), 3 (quad).

use crate::base::error::FEChemError;
use crate::base::mesh::Mesh;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

const CALLER: &str = "Mesh::new";

pub fn read_gmsh_mesh(file_path: &str) -> Result<Mesh, FEChemError> {
    let text = fs::read_to_string(Path::new(file_path)).map_err(|e| FEChemError::MeshFileRead {
        caller: CALLER.to_string(),
        file_path: file_path.to_string(),
        source: e,
    })?;

    let raw = scan_gmsh_sections(&text)?;
    if raw.major_version != 4 {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!(
                "unsupported Gmsh major version {} (need 4.x; got {:?})",
                raw.major_version, raw.version_string
            ),
        });
    }

    let parsed = parse_from_sections_v4(&raw)?;
    assemble_mesh(parsed)
}

struct RawGmshSections<'a> {
    major_version: u32,
    version_string: String,
    entities: Option<Vec<&'a str>>,
    nodes: Vec<&'a str>,
    elements: Vec<&'a str>,
}

fn scan_gmsh_sections(text: &str) -> Result<RawGmshSections<'_>, FEChemError> {
    let mut version_string = String::new();
    let mut entities_block: Option<Vec<&str>> = None;
    let mut nodes_block: Option<Vec<&str>> = None;
    let mut elements_block: Option<Vec<&str>> = None;

    let mut lines = text.lines().map(|l| l.trim()).filter(|l| !l.is_empty());

    while let Some(line) = lines.next() {
        match line {
            "$MeshFormat" => {
                let fmt = lines.next().ok_or_else(|| FEChemError::InvalidGmsh {
                    caller: CALLER.to_string(),
                    message: "unexpected EOF in $MeshFormat".to_string(),
                })?;
                version_string = fmt.split_whitespace().next().unwrap_or("").to_string();
                while let Some(l) = lines.next() {
                    if l == "$EndMeshFormat" {
                        break;
                    }
                }
            }
            "$PhysicalNames" => {
                while let Some(l) = lines.next() {
                    if l == "$EndPhysicalNames" {
                        break;
                    }
                }
            }
            "$Entities" => {
                let mut buf = Vec::new();
                while let Some(l) = lines.next() {
                    if l == "$EndEntities" {
                        break;
                    }
                    buf.push(l);
                }
                entities_block = Some(buf);
            }
            "$Nodes" => {
                let mut buf = Vec::new();
                while let Some(l) = lines.next() {
                    if l == "$EndNodes" {
                        break;
                    }
                    buf.push(l);
                }
                nodes_block = Some(buf);
            }
            "$Elements" => {
                let mut buf = Vec::new();
                while let Some(l) = lines.next() {
                    if l == "$EndElements" {
                        break;
                    }
                    buf.push(l);
                }
                elements_block = Some(buf);
            }
            _ => {}
        }
    }

    let major_version = version_string
        .split('.')
        .next()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);

    let nodes = nodes_block.ok_or_else(|| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: "missing $Nodes section".to_string(),
    })?;
    let elements = elements_block.ok_or_else(|| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: "missing $Elements section".to_string(),
    })?;

    Ok(RawGmshSections {
        major_version,
        version_string,
        entities: entities_block,
        nodes,
        elements,
    })
}

struct ParsedGmsh {
    vert_x: Vec<f64>,
    vert_y: Vec<f64>,
    /// Each 2D cell: ordered CCW vertex indices and raw physical tag from file.
    cells: Vec<(Vec<usize>, i32)>,
    /// Boundary line elements: vertex pair and physical tag.
    lines: Vec<(usize, usize, i32)>,
}

fn parse_from_sections_v4(raw: &RawGmshSections<'_>) -> Result<ParsedGmsh, FEChemError> {
    let entities = raw.entities.as_ref().ok_or_else(|| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: "MSH4 file missing $Entities section (required for physical groups)".to_string(),
    })?;

    let (curve_phys, surface_phys) = parse_entities_v4(entities)?;
    let (vert_x, vert_y, tag_to_vid) = parse_nodes_v4(&raw.nodes)?;
    let (cells, lines) = parse_elements_v4(
        &raw.elements,
        &tag_to_vid,
        &curve_phys,
        &surface_phys,
    )?;

    if cells.is_empty() {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: "no 2D surface elements (entityDim 2, types 2 or 3) found".to_string(),
        });
    }

    Ok(ParsedGmsh {
        vert_x,
        vert_y,
        cells,
        lines,
    })
}

/// `$Entities`: map curve / surface Gmsh entity tags to a representative physical tag (minimum tag if several).
fn parse_entities_v4(block: &[&str]) -> Result<(HashMap<i32, i32>, HashMap<i32, i32>), FEChemError> {
    if block.is_empty() {
        return Ok((HashMap::new(), HashMap::new()));
    }

    let mut it = block[0].split_whitespace();
    let n_points: usize = it
        .next()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Entities: invalid header {:?}", block.get(0)),
        })?;
    let n_curves: usize = it.next().and_then(|s| s.parse().ok()).ok_or_else(|| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: "$Entities: missing curve count".to_string(),
    })?;
    let n_surfaces: usize = it.next().and_then(|s| s.parse().ok()).ok_or_else(|| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: "$Entities: missing surface count".to_string(),
    })?;
    let n_volumes: usize = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);

    let mut line_idx = 1usize;
    for _ in 0..n_points {
        line_idx += 1;
        if line_idx > block.len() {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "$Entities: truncated point block".to_string(),
            });
        }
    }

    let mut curve_phys = HashMap::new();
    for _ in 0..n_curves {
        let line = block.get(line_idx).ok_or_else(|| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: "$Entities: truncated curve block".to_string(),
        })?;
        let (tag, phys) = parse_entity_line_phys_tag(line)?;
        curve_phys.insert(tag, phys);
        line_idx += 1;
    }

    let mut surface_phys = HashMap::new();
    for _ in 0..n_surfaces {
        let line = block.get(line_idx).ok_or_else(|| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: "$Entities: truncated surface block".to_string(),
        })?;
        let (tag, phys) = parse_entity_line_phys_tag(line)?;
        surface_phys.insert(tag, phys);
        line_idx += 1;
    }

    for _ in 0..n_volumes {
        line_idx += 1;
        if line_idx > block.len() {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "$Entities: truncated volume block".to_string(),
            });
        }
    }

    Ok((curve_phys, surface_phys))
}

fn parse_entity_line_phys_tag(line: &str) -> Result<(i32, i32), FEChemError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 8 {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Entities: line too short: {line}"),
        });
    }
    let tag: i32 = tokens[0].parse().map_err(|_| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: format!("$Entities: bad entity tag in {line}"),
    })?;

    let mut idx = 7usize;
    let num_phys: usize = tokens[idx].parse().map_err(|_| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: format!("$Entities: bad numPhysicalTags in {line}"),
    })?;
    idx += 1;

    let phys = if num_phys == 0 {
        0
    } else {
        if tokens.len() < idx + num_phys {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("$Entities: missing physical tags in {line}"),
            });
        }
        let mut pmin: i32 = i32::MAX;
        for _ in 0..num_phys {
            let p: i32 = tokens[idx].parse().map_err(|_| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("$Entities: bad physical tag in {line}"),
            })?;
            pmin = pmin.min(p);
            idx += 1;
        }
        if pmin == i32::MAX {
            0
        } else {
            pmin
        }
    };

    Ok((tag, phys))
}

fn parse_nodes_v4(block: &[&str]) -> Result<(Vec<f64>, Vec<f64>, HashMap<i32, usize>), FEChemError> {
    if block.is_empty() {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: "$Nodes section is empty".to_string(),
        });
    }

    let header: Vec<&str> = block[0].split_whitespace().collect();
    if header.len() < 4 {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Nodes MSH4: bad header {:?}", block[0]),
        });
    }
    let num_entity_blocks: usize = header[0].parse().map_err(|_| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: format!("$Nodes MSH4: bad numEntityBlocks {:?}", header[0]),
    })?;

    let mut tag_to_xy: HashMap<i32, (f64, f64)> = HashMap::new();
    let mut idx = 1usize;

    for _ in 0..num_entity_blocks {
        let h: Vec<&str> = block
            .get(idx)
            .ok_or_else(|| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "$Nodes MSH4: unexpected EOF in entity block header".to_string(),
            })?
            .split_whitespace()
            .collect();
        idx += 1;

        if h.len() < 4 {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("$Nodes MSH4: bad entity block header {:?}", h.join(" ")),
            });
        }

        let parametric: u8 = h[2].parse().map_err(|_| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Nodes MSH4: bad parametric flag {:?}", h[2]),
        })?;
        let num_in_block: usize = h[3].parse().map_err(|_| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Nodes MSH4: bad numNodesInBlock {:?}", h[3]),
        })?;

        if parametric != 0 {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "parametric Gmsh nodes (parametric=1) are not supported".to_string(),
            });
        }

        let mut tags = Vec::with_capacity(num_in_block);
        for _ in 0..num_in_block {
            let tag_line = block.get(idx).ok_or_else(|| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "$Nodes MSH4: EOF reading node tags".to_string(),
            })?;
            let node_tag: i32 = tag_line.trim().parse().map_err(|_| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("$Nodes MSH4: bad node tag line {tag_line}"),
            })?;
            tags.push(node_tag);
            idx += 1;
        }

        for &node_tag in &tags {
            let coord_line = block.get(idx).ok_or_else(|| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "$Nodes MSH4: EOF reading node coordinates".to_string(),
            })?;
            idx += 1;
            let mut c = coord_line.split_whitespace();
            let x: f64 = c.next().and_then(|s| s.parse().ok()).ok_or_else(|| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("$Nodes MSH4: bad x in {coord_line}"),
            })?;
            let y: f64 = c.next().and_then(|s| s.parse().ok()).ok_or_else(|| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("$Nodes MSH4: bad y in {coord_line}"),
            })?;
            tag_to_xy.insert(node_tag, (x, y));
        }
    }

    if idx != block.len() {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!(
                "$Nodes MSH4: expected {} lines in section, have {}",
                idx,
                block.len()
            ),
        });
    }

    let mut pairs: Vec<(i32, f64, f64)> = tag_to_xy.into_iter().map(|(t, (x, y))| (t, x, y)).collect();
    pairs.sort_by_key(|p| p.0);

    let mut tag_to_vid = HashMap::new();
    let mut vert_x = Vec::with_capacity(pairs.len());
    let mut vert_y = Vec::with_capacity(pairs.len());
    for (vid, (tag, x, y)) in pairs.into_iter().enumerate() {
        tag_to_vid.insert(tag, vid);
        vert_x.push(x);
        vert_y.push(y);
    }

    Ok((vert_x, vert_y, tag_to_vid))
}

fn msh4_element_line_node_count(elm_type: usize) -> Result<usize, FEChemError> {
    match elm_type {
        1 => Ok(2),
        2 => Ok(3),
        3 => Ok(4),
        15 => Ok(1),
        _ => Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!(
                "MSH4 element type {elm_type} is not supported (only linear types 1, 2, 3 and point 15 for skipping)"
            ),
        }),
    }
}

fn parse_elements_v4(
    block: &[&str],
    tag_to_vid: &HashMap<i32, usize>,
    curve_phys: &HashMap<i32, i32>,
    surface_phys: &HashMap<i32, i32>,
) -> Result<(Vec<(Vec<usize>, i32)>, Vec<(usize, usize, i32)>), FEChemError> {
    if block.is_empty() {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: "$Elements section is empty".to_string(),
        });
    }

    let header: Vec<&str> = block[0].split_whitespace().collect();
    if header.len() < 4 {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Elements MSH4: bad header {:?}", block[0]),
        });
    }
    let num_entity_blocks: usize = header[0].parse().map_err(|_| FEChemError::InvalidGmsh {
        caller: CALLER.to_string(),
        message: format!("$Elements MSH4: bad numEntityBlocks {:?}", header[0]),
    })?;

    let mut cells = Vec::new();
    let mut lines = Vec::new();
    let mut idx = 1usize;

    for _ in 0..num_entity_blocks {
        let h: Vec<&str> = block
            .get(idx)
            .ok_or_else(|| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "$Elements MSH4: EOF block header".to_string(),
            })?
            .split_whitespace()
            .collect();
        idx += 1;

        if h.len() < 4 {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("$Elements MSH4: bad block header {:?}", h.join(" ")),
            });
        }

        let entity_dim: usize = h[0].parse().map_err(|_| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Elements MSH4: bad entityDim {:?}", h[0]),
        })?;
        let entity_tag: i32 = h[1].parse().map_err(|_| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Elements MSH4: bad entityTag {:?}", h[1]),
        })?;
        let elm_type: usize = h[2].parse().map_err(|_| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Elements MSH4: bad elementType {:?}", h[2]),
        })?;
        let num_elm: usize = h[3].parse().map_err(|_| FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!("$Elements MSH4: bad numElementsInBlock {:?}", h[3]),
        })?;

        let nn = match msh4_element_line_node_count(elm_type) {
            Ok(n) => n,
            Err(_) => {
                return Err(FEChemError::InvalidGmsh {
                    caller: CALLER.to_string(),
                    message: format!(
                        "MSH4: unsupported element type {elm_type} in entityDim {entity_dim} block (use linear 1/2/3 only)"
                    ),
                });
            }
        };

        for _ in 0..num_elm {
            let line = block.get(idx).ok_or_else(|| FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: "$Elements MSH4: EOF element line".to_string(),
            })?;
            idx += 1;

            let nums: Vec<i32> = line.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            if nums.len() < 1 + nn {
                return Err(FEChemError::InvalidGmsh {
                    caller: CALLER.to_string(),
                    message: format!("$Elements MSH4: short line {line}"),
                });
            }

            let node_start = 1;
            match (entity_dim, elm_type) {
                (1, 1) => {
                    let t0 = nums[node_start];
                    let t1 = nums[node_start + 1];
                    let v0 = *tag_to_vid.get(&t0).ok_or_else(|| FEChemError::InvalidGmsh {
                        caller: CALLER.to_string(),
                        message: format!("line element references unknown node tag {t0}"),
                    })?;
                    let v1 = *tag_to_vid.get(&t1).ok_or_else(|| FEChemError::InvalidGmsh {
                        caller: CALLER.to_string(),
                        message: format!("line element references unknown node tag {t1}"),
                    })?;
                    let phys = curve_phys.get(&entity_tag).copied().unwrap_or(0);
                    lines.push((v0, v1, phys));
                }
                (2, 2) => {
                    let mut v = Vec::with_capacity(3);
                    for k in 0..3 {
                        let t = nums[node_start + k];
                        v.push(*tag_to_vid.get(&t).ok_or_else(|| FEChemError::InvalidGmsh {
                            caller: CALLER.to_string(),
                            message: format!("triangle references unknown node tag {t}"),
                        })?);
                    }
                    let phys = surface_phys.get(&entity_tag).copied().unwrap_or(0);
                    cells.push((v, phys));
                }
                (2, 3) => {
                    let mut v = Vec::with_capacity(4);
                    for k in 0..4 {
                        let t = nums[node_start + k];
                        v.push(*tag_to_vid.get(&t).ok_or_else(|| FEChemError::InvalidGmsh {
                            caller: CALLER.to_string(),
                            message: format!("quad references unknown node tag {t}"),
                        })?);
                    }
                    let phys = surface_phys.get(&entity_tag).copied().unwrap_or(0);
                    cells.push((v, phys));
                }
                (_, 15) => {}
                _ => {}
            }
        }
    }

    if idx != block.len() {
        return Err(FEChemError::InvalidGmsh {
            caller: CALLER.to_string(),
            message: format!(
                "$Elements MSH4: section line count mismatch (parsed to {}, len {})",
                idx,
                block.len()
            ),
        });
    }

    Ok((cells, lines))
}

fn polygon_signed_area(indices: &[usize], vx: &[f64], vy: &[f64]) -> f64 {
    let n = indices.len();
    let mut a = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let xi = vx[indices[i]];
        let yi = vy[indices[i]];
        let xj = vx[indices[j]];
        let yj = vy[indices[j]];
        a += xi * yj - xj * yi;
    }
    0.5 * a
}

fn order_cell_ccw(v: &mut Vec<usize>, vx: &[f64], vy: &[f64]) {
    if v.len() < 3 {
        return;
    }
    if polygon_signed_area(v, vx, vy) < 0.0 {
        v.reverse();
    }
}

/// Map raw physical tags to consecutive region indices `0..n-1`.
fn consecutive_phys_remap(tags: &[i32]) -> HashMap<i32, usize> {
    let mut uniq: Vec<i32> = tags.iter().copied().collect::<HashSet<_>>().into_iter().collect();
    uniq.sort();
    uniq.iter().enumerate().map(|(i, &t)| (t, i)).collect()
}

fn assemble_mesh(parsed: ParsedGmsh) -> Result<Mesh, FEChemError> {
    let ParsedGmsh {
        vert_x,
        vert_y,
        mut cells,
        lines,
    } = parsed;

    for (v, _) in cells.iter_mut() {
        order_cell_ccw(v, &vert_x, &vert_y);
    }

    let num_elm2d = cells.len();
    let mut elm2d_node = Vec::with_capacity(num_elm2d);
    let mut elm2d_node_id = Vec::with_capacity(num_elm2d);
    let mut raw_reg2d_tags = Vec::with_capacity(num_elm2d);

    for (ci, (verts, phys)) in cells.iter().enumerate() {
        let area = polygon_signed_area(verts, &vert_x, &vert_y).abs();
        if area <= 1e-30 {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("degenerate 2D element at index {ci}"),
            });
        }
        elm2d_node.push(verts.len());
        elm2d_node_id.push(verts.clone());
        raw_reg2d_tags.push(*phys);
    }

    let reg2d_remap = consecutive_phys_remap(&raw_reg2d_tags);
    let num_reg2d = reg2d_remap.len();
    let mut reg2d_elem_id: Vec<Vec<usize>> = vec![Vec::new(); num_reg2d];
    for (ei, &phys) in raw_reg2d_tags.iter().enumerate() {
        reg2d_elem_id[reg2d_remap[&phys]].push(ei);
    }

    let num_elm1d = lines.len();
    let mut elm1d_node = Vec::with_capacity(num_elm1d);
    let mut elm1d_node_id = Vec::with_capacity(num_elm1d);
    let mut raw_reg1d_tags = Vec::with_capacity(num_elm1d);

    for (ei, &(v0, v1, phys)) in lines.iter().enumerate() {
        let dx = vert_x[v1] - vert_x[v0];
        let dy = vert_y[v1] - vert_y[v0];
        if dx * dx + dy * dy <= 1e-30 {
            return Err(FEChemError::InvalidGmsh {
                caller: CALLER.to_string(),
                message: format!("degenerate 1D element at index {ei}"),
            });
        }
        elm1d_node.push(2);
        elm1d_node_id.push(vec![v0, v1]);
        raw_reg1d_tags.push(phys);
    }

    let reg1d_remap = consecutive_phys_remap(&raw_reg1d_tags);
    let num_reg1d = reg1d_remap.len();
    let mut reg1d_elem_id: Vec<Vec<usize>> = vec![Vec::new(); num_reg1d];
    for (ei, &phys) in raw_reg1d_tags.iter().enumerate() {
        reg1d_elem_id[reg1d_remap[&phys]].push(ei);
    }

    Ok(Mesh {
        num_node: vert_x.len(),
        node_x: vert_x,
        node_y: vert_y,
        num_elm2d,
        num_reg2d,
        elm2d_node_num: elm2d_node,
        elm2d_node_id,
        reg2d_elem_id,
        num_elm1d,
        num_reg1d,
        elm1d_node_num: elm1d_node,
        elm1d_node_id,
        reg1d_elem_id,
    })
}

#[cfg(test)]
mod tests {
    use super::{consecutive_phys_remap, read_gmsh_mesh};
    use crate::base::geom_bnd::Boundary;
    use crate::base::geom_dom::Domain;
    use crate::base::mesh::Mesh;

    fn assert_mesh_topology(mesh: &Mesh) {
        assert_eq!(mesh.node_x.len(), mesh.num_node);
        assert_eq!(mesh.node_y.len(), mesh.num_node);
        assert_eq!(mesh.elm2d_node_num.len(), mesh.num_elm2d);
        assert_eq!(mesh.elm2d_node_id.len(), mesh.num_elm2d);
        assert_eq!(mesh.elm1d_node_num.len(), mesh.num_elm1d);
        assert_eq!(mesh.elm1d_node_id.len(), mesh.num_elm1d);
        assert_eq!(mesh.reg2d_elem_id.len(), mesh.num_reg2d);
        assert_eq!(mesh.reg1d_elem_id.len(), mesh.num_reg1d);

        let reg2d_count: usize = mesh.reg2d_elem_id.iter().map(|r| r.len()).sum();
        assert_eq!(reg2d_count, mesh.num_elm2d);

        let reg1d_count: usize = mesh.reg1d_elem_id.iter().map(|r| r.len()).sum();
        assert_eq!(reg1d_count, mesh.num_elm1d);

        for (ei, nodes) in mesh.elm2d_node_id.iter().enumerate() {
            assert_eq!(nodes.len(), mesh.elm2d_node_num[ei]);
            assert!(nodes.len() == 3 || nodes.len() == 4);
            for &nid in nodes {
                assert!(nid < mesh.num_node);
            }
        }

        for (ei, nodes) in mesh.elm1d_node_id.iter().enumerate() {
            assert_eq!(nodes.len(), 2);
            assert_eq!(mesh.elm1d_node_num[ei], 2);
            for &nid in nodes {
                assert!(nid < mesh.num_node);
            }
        }
    }

    #[test]
    fn consecutive_phys_remap_non_contiguous() {
        let tags = vec![10, 20, 10, 30, 20];
        let remap = consecutive_phys_remap(&tags);
        assert_eq!(remap.len(), 3);
        assert_eq!(remap[&10], 0);
        assert_eq!(remap[&20], 1);
        assert_eq!(remap[&30], 2);
    }

    #[test]
    fn gmsh_41_square_tri_mesh() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/square_tri.msh");
        let mesh = read_gmsh_mesh(path).unwrap();

        assert!(mesh.num_node > 4);
        assert!(mesh.num_elm2d > 0);
        assert_eq!(mesh.num_reg2d, 1);
        assert_eq!(mesh.num_reg1d, 4);
        assert_mesh_topology(&mesh);

        let dom = Domain::new(0, &mesh, 0).unwrap();
        assert_eq!(dom.num_elem, mesh.num_elm2d);
        assert!(dom.num_node > 0);

        for reg_id in 0..mesh.num_reg1d {
            let _bnd = Boundary::new(reg_id, &mesh, &dom, reg_id).unwrap();
        }
    }
}
