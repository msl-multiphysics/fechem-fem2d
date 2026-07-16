import gmsh

# Square with bottom-left and top-right corner rectangles
#
# Boundaries and Interfaces
# 1 - bottom-left outer (left + bottom)
# 2 - middle outer (bottom + right)
# 3 - top-right outer (right + top)
# 4 - middle outer (top + left)
# 5 - bottom-left interface
# 6 - top-right interface
#
# Domains
# 1 - middle
# 2 - bottom-left
# 3 - top-right
#

# initialize gmsh
gmsh.initialize()
gmsh.model.add("gmsh_threereg")

# create points
size = 0.05
p1 = gmsh.model.geo.addPoint(0.00, 0.00, 0, size, 1)
p2 = gmsh.model.geo.addPoint(0.75, 0.00, 0, size, 2)
p3 = gmsh.model.geo.addPoint(1.00, 0.00, 0, size, 3)
p4 = gmsh.model.geo.addPoint(1.00, 0.75, 0, size, 4)
p5 = gmsh.model.geo.addPoint(1.00, 1.00, 0, size, 5)
p6 = gmsh.model.geo.addPoint(0.25, 1.00, 0, size, 6)
p7 = gmsh.model.geo.addPoint(0.00, 1.00, 0, size, 7)
p8 = gmsh.model.geo.addPoint(0.00, 0.25, 0, size, 8)
p9 = gmsh.model.geo.addPoint(0.75, 0.25, 0, size, 9)
p10 = gmsh.model.geo.addPoint(0.25, 0.75, 0, size, 10)

# create curves
c1 = gmsh.model.geo.addLine(p1, p2, 1)    # bottom-left bottom
c2 = gmsh.model.geo.addLine(p2, p3, 2)    # middle bottom
c3 = gmsh.model.geo.addLine(p3, p4, 3)    # middle right
c4 = gmsh.model.geo.addLine(p4, p5, 4)    # top-right right
c5 = gmsh.model.geo.addLine(p5, p6, 5)    # top-right top
c6 = gmsh.model.geo.addLine(p6, p7, 6)    # middle top
c7 = gmsh.model.geo.addLine(p7, p8, 7)    # middle left
c8 = gmsh.model.geo.addLine(p8, p1, 8)    # bottom-left left
c9 = gmsh.model.geo.addLine(p2, p9, 9)    # bottom-left interface (vertical)
c10 = gmsh.model.geo.addLine(p9, p8, 10)  # bottom-left interface (horizontal)
c11 = gmsh.model.geo.addLine(p6, p10, 11) # top-right interface (vertical)
c12 = gmsh.model.geo.addLine(p10, p4, 12) # top-right interface (horizontal)

# create loops
l1 = gmsh.model.geo.addCurveLoop([c2, c3, -c12, -c11, c6, c7, -c10, -c9], 1)  # middle
l2 = gmsh.model.geo.addCurveLoop([c1, c9, c10, c8], 2)                        # bottom-left
l3 = gmsh.model.geo.addCurveLoop([c4, c5, c11, c12], 3)                       # top-right

# create surfaces
s1 = gmsh.model.geo.addPlaneSurface([l1], 1)
s2 = gmsh.model.geo.addPlaneSurface([l2], 2)
s3 = gmsh.model.geo.addPlaneSurface([l3], 3)

# synchronize
gmsh.model.geo.synchronize()

# add physical groups
# FEChem will use physical groups to number regions and boundaries
# for a single domain, the FEChem indices are the gmsh indices minus 1
gmsh.model.addPhysicalGroup(1, [c8, c1], 1)   # bottom-left outer     - 1
gmsh.model.addPhysicalGroup(1, [c2, c3], 2)   # middle outer          - 2
gmsh.model.addPhysicalGroup(1, [c4, c5], 3)   # top-right outer       - 3
gmsh.model.addPhysicalGroup(1, [c6, c7], 4)   # middle outer          - 4
gmsh.model.addPhysicalGroup(1, [c9, c10], 5)  # bottom-left interface - 5
gmsh.model.addPhysicalGroup(1, [c11, c12], 6) # top-right interface   - 6
gmsh.model.addPhysicalGroup(2, [s1], 1)       # middle                - 1
gmsh.model.addPhysicalGroup(2, [s2], 2)       # bottom-left           - 2
gmsh.model.addPhysicalGroup(2, [s3], 3)       # top-right             - 3

# generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("gmsh_threereg.msh")
gmsh.finalize()
