import gmsh

# Channel with a circular cylinder for a von Karman vortex street
# Has smaller tri elements near the cylinder and in its downstream wake
#
# Boundaries
# 1 - left
# 2 - right
# 3 - bottom
# 4 - top
# 5 - cylinder
#
# Domains
# 1 - domain
#

# initialize gmsh
gmsh.initialize()
gmsh.model.add("gmsh_vortex")

# mesh parameters
length = 2.2  # channel length
height = 0.41  # channel height
cylinder_x = 0.2  # cylinder center x-coordinate
cylinder_y = 0.2  # cylinder center y-coordinate
cylinder_radius = 0.05  # cylinder radius
size_cylinder = 0.01  # size of elements near the cylinder
size_wake = 0.03  # size of elements in the downstream wake
size_far = 0.10  # size of elements far from the cylinder and wake
ref_start = 0.05  # distance from cylinder to start coarsening
ref_end = 0.3  # distance from cylinder to end coarsening

# create points
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0, size_far, 1)
p2 = gmsh.model.geo.addPoint(length, 0.0, 0, size_far, 2)
p3 = gmsh.model.geo.addPoint(length, height, 0, size_far, 3)
p4 = gmsh.model.geo.addPoint(0.0, height, 0, size_far, 4)
p5 = gmsh.model.geo.addPoint(cylinder_x, cylinder_y, 0, size_cylinder, 5)
p6 = gmsh.model.geo.addPoint(
    cylinder_x + cylinder_radius, cylinder_y, 0, size_cylinder, 6
)
p7 = gmsh.model.geo.addPoint(
    cylinder_x, cylinder_y + cylinder_radius, 0, size_cylinder, 7
)
p8 = gmsh.model.geo.addPoint(
    cylinder_x - cylinder_radius, cylinder_y, 0, size_cylinder, 8
)
p9 = gmsh.model.geo.addPoint(
    cylinder_x, cylinder_y - cylinder_radius, 0, size_cylinder, 9
)

# create curves
c1 = gmsh.model.geo.addLine(p4, p1, 1)  # left   - 1
c2 = gmsh.model.geo.addLine(p2, p3, 2)  # right  - 2
c3 = gmsh.model.geo.addLine(p1, p2, 3)  # bottom - 3
c4 = gmsh.model.geo.addLine(p3, p4, 4)  # top    - 4
c5 = gmsh.model.geo.addCircleArc(p6, p5, p7, 5)  # cylinder - 5
c6 = gmsh.model.geo.addCircleArc(p7, p5, p8, 6)  # cylinder - 6
c7 = gmsh.model.geo.addCircleArc(p8, p5, p9, 7)  # cylinder - 7
c8 = gmsh.model.geo.addCircleArc(p9, p5, p6, 8)  # cylinder - 8

# create loops
l1 = gmsh.model.geo.addCurveLoop([c1, c3, c2, c4], 1)
l2 = gmsh.model.geo.addCurveLoop([c5, c6, c7, c8], 2)

# create surface
s1 = gmsh.model.geo.addPlaneSurface([l1, l2], 1)

# synchronize
gmsh.model.geo.synchronize()

# add physical groups
# FEChem will use physical groups to number regions and boundaries
# for a single domain, the FEChem indices are the gmsh indices minus 1
gmsh.model.addPhysicalGroup(1, [c1], 1)  # left   - 1
gmsh.model.addPhysicalGroup(1, [c2], 2)  # right  - 2
gmsh.model.addPhysicalGroup(1, [c3], 3)  # bottom - 3
gmsh.model.addPhysicalGroup(1, [c4], 4)  # top    - 4
gmsh.model.addPhysicalGroup(1, [c5, c6, c7, c8], 5)  # cylinder - 5
gmsh.model.addPhysicalGroup(2, [s1], 1)  # domain - 1

# cylinder-refined mesh-size field
# prescribe element size from distance to the cylinder:
#   distance <= ref_start  -> size_cylinder
#   distance >= ref_end    -> size_far
#   in between             -> smooth transition

# distance from each mesh location to the cylinder boundary
distance_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance_field, "CurvesList", [c5, c6, c7, c8])
gmsh.model.mesh.field.setNumber(distance_field, "Sampling", 100)  # boundary sample points

# map distance field to target element size
threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", size_cylinder)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", size_far)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", ref_start)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", ref_end)

# refine the wake downstream of the cylinder
wake_field = gmsh.model.mesh.field.add("Box")
gmsh.model.mesh.field.setNumber(wake_field, "VIn", size_wake)
gmsh.model.mesh.field.setNumber(wake_field, "VOut", size_far)
gmsh.model.mesh.field.setNumber(wake_field, "XMin", cylinder_x)
gmsh.model.mesh.field.setNumber(wake_field, "XMax", length)
gmsh.model.mesh.field.setNumber(
    wake_field, "YMin", cylinder_y - 2.5 * cylinder_radius
)
gmsh.model.mesh.field.setNumber(
    wake_field, "YMax", cylinder_y + 2.5 * cylinder_radius
)
gmsh.model.mesh.field.setNumber(wake_field, "Thickness", 0.08)

# use the finest size prescribed by the cylinder and wake fields
minimum_field = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(
    minimum_field, "FieldsList", [threshold_field, wake_field]
)
gmsh.model.mesh.field.setAsBackgroundMesh(minimum_field)

# let the background field control sizing, not point/curve defaults
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

# generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("gmsh_vortex.msh")
gmsh.finalize()
