import gmsh

# Square with uniformly sized tri elements divided into three triangles
# Has a larger diagonal going from the lower left to the upper right
# Has a smaller diagonal going from the lower right to the center
#
# Boundaries and Interfaces
# 1 - left
# 2 - right
# 3 - bottom
# 4 - top
# 5 - large diagonal (lower-left half)
# 6 - large diagonal (upper-right half)
# 7 - small diagonal
#
# Domains
# 1 - upper-left triangle
# 2 - right triangle
# 3 - bottom triangle
#

# initialize gmsh
gmsh.initialize()
gmsh.model.add("gmsh_uniform_multi")

# create points
size = 0.05
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0, size, 1)
p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0, size, 2)
p3 = gmsh.model.geo.addPoint(1.0, 1.0, 0, size, 3)
p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0, size, 4)
p5 = gmsh.model.geo.addPoint(0.5, 0.5, 0, size, 5)  # square center

# create curves
c1 = gmsh.model.geo.addLine(p4, p1, 1)  # left           - 1
c2 = gmsh.model.geo.addLine(p2, p3, 2)  # right          - 2
c3 = gmsh.model.geo.addLine(p1, p2, 3)  # bottom         - 3
c4 = gmsh.model.geo.addLine(p3, p4, 4)  # top            - 4
c5 = gmsh.model.geo.addLine(p1, p5, 5)  # large diagonal - 5 (lower-left half)
c6 = gmsh.model.geo.addLine(p5, p3, 6)  # large diagonal - 6 (upper-right half)
c7 = gmsh.model.geo.addLine(p2, p5, 7)  # small diagonal - 7

# create loops
l1 = gmsh.model.geo.addCurveLoop([c1, c5, c6, c4], 1)  # upper-left triangle
l2 = gmsh.model.geo.addCurveLoop([c2, -c6, -c7], 2)    # right triangle
l3 = gmsh.model.geo.addCurveLoop([c3, c7, -c5], 3)     # bottom triangle

# create surfaces
s1 = gmsh.model.geo.addPlaneSurface([l1], 1)
s2 = gmsh.model.geo.addPlaneSurface([l2], 2)
s3 = gmsh.model.geo.addPlaneSurface([l3], 3)

# synchronize
gmsh.model.geo.synchronize()

# add physical groups
# FEChem will use physical groups to number regions and boundaries
# for a single domain, the FEChem indices are the gmsh indices minus 1
gmsh.model.addPhysicalGroup(1, [c1], 1)  # left                              - 1
gmsh.model.addPhysicalGroup(1, [c2], 2)  # right                             - 2
gmsh.model.addPhysicalGroup(1, [c3], 3)  # bottom                            - 3
gmsh.model.addPhysicalGroup(1, [c4], 4)  # top                               - 4
gmsh.model.addPhysicalGroup(1, [c5], 5)  # large diagonal (lower-left half)  - 5
gmsh.model.addPhysicalGroup(1, [c6], 6)  # large diagonal (upper-right half) - 6
gmsh.model.addPhysicalGroup(1, [c7], 7)  # small diagonal                    - 7
gmsh.model.addPhysicalGroup(2, [s1], 1)  # upper-left triangle               - 1
gmsh.model.addPhysicalGroup(2, [s2], 2)  # right triangle                    - 2
gmsh.model.addPhysicalGroup(2, [s3], 3)  # bottom triangle                   - 3

# generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("gmsh_uniform_multi.msh")
gmsh.finalize()
