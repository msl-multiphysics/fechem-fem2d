import gmsh

# Channel with bump surrounded by solid regions.
#
# Boundaries and Interfaces
# 1 - inlet
# 2 - outlet
# 3 - bottom solid
# 4 - top solid
# 5 - bottom solid-channel interface
# 6 - top solid-channel interface
#
# Domains
# 1 - channel
# 2 - bottom solid
# 3 - top solid
#

# initialize gmsh
gmsh.initialize()
gmsh.model.add("gmsh_heater")

# create points
size_solid = 0.05
size_fluid = 0.02
p1  = gmsh.model.geo.addPoint(0.00, 0.00, 0, size_fluid, 1)
p2  = gmsh.model.geo.addPoint(1.50, 0.00, 0, size_fluid, 2)
p3  = gmsh.model.geo.addPoint(1.50, 0.50, 0, size_fluid, 3)
p4  = gmsh.model.geo.addPoint(2.00, 0.50, 0, size_fluid, 4)
p5  = gmsh.model.geo.addPoint(2.00, 0.00, 0, size_fluid, 5)
p6  = gmsh.model.geo.addPoint(3.50, 0.00, 0, size_fluid, 6)

p7  = gmsh.model.geo.addPoint(0.00, 0.50, 0, size_fluid, 7)
p8  = gmsh.model.geo.addPoint(1.00, 0.50, 0, size_fluid, 8)
p9  = gmsh.model.geo.addPoint(1.00, 1.00, 0, size_fluid, 9)
p10 = gmsh.model.geo.addPoint(2.50, 1.00, 0, size_fluid, 10)
p11 = gmsh.model.geo.addPoint(2.50, 0.50, 0, size_fluid, 11)
p12 = gmsh.model.geo.addPoint(3.50, 0.50, 0, size_fluid, 12)

p13 = gmsh.model.geo.addPoint(0.00, -0.50, 0, size_solid, 13)
p14 = gmsh.model.geo.addPoint(3.50, -0.50, 0, size_solid, 14)

p15 = gmsh.model.geo.addPoint(0.00, 1.50, 0, size_solid, 15)
p16 = gmsh.model.geo.addPoint(3.50, 1.50, 0, size_solid, 16)

# create curves
c1 = gmsh.model.geo.addLine(p1, p2, 1)
c2 = gmsh.model.geo.addLine(p2, p3, 2)
c3 = gmsh.model.geo.addLine(p3, p4, 3)
c4 = gmsh.model.geo.addLine(p4, p5, 4)
c5 = gmsh.model.geo.addLine(p5, p6, 5)

c6  = gmsh.model.geo.addLine( p7,  p8, 6)
c7  = gmsh.model.geo.addLine( p8,  p9, 7)
c8  = gmsh.model.geo.addLine( p9, p10, 8)
c9  = gmsh.model.geo.addLine(p10, p11, 9)
c10 = gmsh.model.geo.addLine(p11, p12, 10)

c11 = gmsh.model.geo.addLine(p1,  p7, 11)
c12 = gmsh.model.geo.addLine(p6, p12, 12)

c13 = gmsh.model.geo.addLine( p1, p13, 13)
c14 = gmsh.model.geo.addLine(p13, p14, 14)
c15 = gmsh.model.geo.addLine(p14,  p6, 15)

c16 = gmsh.model.geo.addLine( p7, p15, 16)
c17 = gmsh.model.geo.addLine(p15, p16, 17)
c18 = gmsh.model.geo.addLine(p16, p12, 18)

# create loops
l1 = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4, c5, c12, -c10, -c9, -c8, -c7, -c6, -c11], 1)  # channel
l2 = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4, c5, -c15, -c14, -c13], 2)  # bottom solid
l3 = gmsh.model.geo.addCurveLoop([c6, c7, c8, c9, c10, -c18, -c17, -c16], 3)  # top solid

# create surfaces
s1 = gmsh.model.geo.addPlaneSurface([l1], 1)
s2 = gmsh.model.geo.addPlaneSurface([l2], 2)
s3 = gmsh.model.geo.addPlaneSurface([l3], 3)

# synchronize
gmsh.model.geo.synchronize()

# add physical groups
# FEChem will use physical groups to number regions and boundaries
# for a single domain, the FEChem indices are the gmsh indices minus 1
gmsh.model.addPhysicalGroup(1, [c11], 1)  # inlet
gmsh.model.addPhysicalGroup(1, [c12], 2)  # outlet
gmsh.model.addPhysicalGroup(1, [c13, c14, c15], 3)  # bottom solid-channel interface
gmsh.model.addPhysicalGroup(1, [c16, c17, c18], 4)  # top solid-channel interface
gmsh.model.addPhysicalGroup(1, [c1, c2, c3, c4, c5], 5)  # bottom solid-channel
gmsh.model.addPhysicalGroup(1, [c6, c7, c8, c9, c10], 6)  # top solid-channel
gmsh.model.addPhysicalGroup(2, [s1], 1)  # channel
gmsh.model.addPhysicalGroup(2, [s2], 2)  # bottom solid
gmsh.model.addPhysicalGroup(2, [s3], 3)  # top solid

# generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("gmsh_heater.msh")
gmsh.finalize()
