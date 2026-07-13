import gmsh

# Square with uniformly sized tri elements
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
gmsh.model.add("gmsh_uniform")

# create points
size = 0.05
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0, size, 1)
p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0, size, 2)
p3 = gmsh.model.geo.addPoint(1.0, 1.0, 0, size, 3)
p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0, size, 4)

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

# generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("gmsh_uniform.msh")
gmsh.finalize()
