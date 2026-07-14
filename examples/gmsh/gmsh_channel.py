import math

import gmsh

# Z-shaped tube rotated 45 degrees counterclockwise
#
# Boundaries
# 1 - inlet
# 2 - outlet
# 3 - left wall
# 4 - right wall
#
# Domains
# 1 - domain
#

def offset_intersection(point, direction_before, direction_after, distance):
    """Intersect two lines offset by distance to the left of their directions."""
    x, y = point
    ux1, uy1 = direction_before
    ux2, uy2 = direction_after
    x1 = x - distance * uy1
    y1 = y + distance * ux1
    x2 = x - distance * uy2
    y2 = y + distance * ux2
    cross = ux1 * uy2 - uy1 * ux2
    t = ((x2 - x1) * uy2 - (y2 - y1) * ux2) / cross
    return x1 + t * ux1, y1 + t * uy1


# initialize gmsh
gmsh.initialize()
gmsh.model.add("gmsh_channel")

# geometry and mesh parameters
tube_width = 0.2
segment_length = 0.5
mesh_size = 0.01
diagonal = segment_length / math.sqrt(2.0)

# The unrotated centerline is a Z. After a 45-degree counterclockwise
# rotation, its first and last sections rise diagonally and its middle
# section points straight down.
centerline = [
    (0.0, 0.0),
    (diagonal, diagonal),
    (diagonal, diagonal - segment_length),
    (2.0 * diagonal, 2.0 * diagonal - segment_length),
]

directions = []
for start, end in zip(centerline[:-1], centerline[1:]):
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = math.hypot(dx, dy)
    directions.append((dx / length, dy / length))

half_width = tube_width / 2.0
left = [
    (
        centerline[0][0] - half_width * directions[0][1],
        centerline[0][1] + half_width * directions[0][0],
    )
]
right = [
    (
        centerline[0][0] + half_width * directions[0][1],
        centerline[0][1] - half_width * directions[0][0],
    )
]

for i in range(1, len(centerline) - 1):
    left.append(
        offset_intersection(
            centerline[i], directions[i - 1], directions[i], half_width
        )
    )
    right.append(
        offset_intersection(
            centerline[i], directions[i - 1], directions[i], -half_width
        )
    )

left.append(
    (
        centerline[-1][0] - half_width * directions[-1][1],
        centerline[-1][1] + half_width * directions[-1][0],
    )
)
right.append(
    (
        centerline[-1][0] + half_width * directions[-1][1],
        centerline[-1][1] - half_width * directions[-1][0],
    )
)

# create the tube perimeter
perimeter = left + list(reversed(right))
points = [
    gmsh.model.geo.addPoint(x, y, 0.0, mesh_size)
    for x, y in perimeter
]
curves = [
    gmsh.model.geo.addLine(points[i], points[(i + 1) % len(points)])
    for i in range(len(points))
]

loop = gmsh.model.geo.addCurveLoop(curves)
surface = gmsh.model.geo.addPlaneSurface([loop])

# synchronize
gmsh.model.geo.synchronize()

# add physical groups
# FEChem will use physical groups to number regions and boundaries
# for a single domain, the FEChem indices are the gmsh indices minus 1
gmsh.model.addPhysicalGroup(1, [curves[-1]], 1)  # inlet
gmsh.model.addPhysicalGroup(1, [curves[3]], 2)  # outlet
gmsh.model.addPhysicalGroup(1, curves[:3], 3)  # left wall
gmsh.model.addPhysicalGroup(1, curves[4:7], 4)  # right wall
gmsh.model.addPhysicalGroup(2, [surface], 1)  # domain

# generate mesh
gmsh.model.mesh.generate(2)
gmsh.write("gmsh_channel.msh")
gmsh.finalize()
