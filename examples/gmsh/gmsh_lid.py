import gmsh

# Square with smaller tri elements near walls
#
# Boundaries
# 1 - left
# 2 - right
# 3 - bottom
# 4 - top
#
# Domains 
# 1 - domain
#

# initialize gmsh
gmsh.initialize()
gmsh.model.add("gmsh_lid")

# mesh parameters
size_wall = 0.02  # size of elements near walls
size_center = 0.05  # size of elements in the center
ref_start = 0.04  # distance from wall to start refining
ref_end = 0.25  # distance from wall to end refining

# create points
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0, size_center, 1)
p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0, size_center, 2)
p3 = gmsh.model.geo.addPoint(1.0, 1.0, 0, size_center, 3)
p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0, size_center, 4)

# create curves
c1 = gmsh.model.geo.addLine(p4, p1, 1)  # left   - 1
c2 = gmsh.model.geo.addLine(p2, p3, 2)  # right  - 2
c3 = gmsh.model.geo.addLine(p1, p2, 3)  # bottom - 3
c4 = gmsh.model.geo.addLine(p3, p4, 4)  # top    - 4

# create loop
l1 = gmsh.model.geo.addCurveLoop([c1, c3, c2, c4], 1)

# create surface
s1 = gmsh.model.geo.addPlaneSurface([l1], 1)

# synchronize
gmsh.model.geo.synchronize()

# add physical groups
# FEChem will use physical groups to number regions and boundaries
# for a single domain, the FEChem indices are the gmsh indices minus 1
gmsh.model.addPhysicalGroup(1, [c1], 1)  # left   - 1
gmsh.model.addPhysicalGroup(1, [c2], 2)  # right  - 2
gmsh.model.addPhysicalGroup(1, [c3], 3)  # bottom - 3
gmsh.model.addPhysicalGroup(1, [c4], 4)  # top    - 4
gmsh.model.addPhysicalGroup(2, [s1], 1)  # domain - 1

# wall-refined mesh-size field
# prescribe element size from distance to the cavity walls:
#   distance <= ref_start  -> size_wall
#   distance >= ref_end    -> size_center
#   in between             -> smooth transition

# distance from each mesh location to the four boundary curves
distance_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance_field, "CurvesList", [c1, c2, c3, c4])
gmsh.model.mesh.field.setNumber(distance_field, "Sampling", 200)  # boundary sample points

# map distance field to target element size
threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", size_wall)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", size_center)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", ref_start)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", ref_end)
gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

# let the background field control sizing, not point/curve defaults
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

# generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("gmsh_lid.msh")
gmsh.finalize()
