 _____ _____ ____ _
|  ___| ____/ ___| |__   ___ _ __ ___
| |_  |  _|| |   | '_ \ / _ \ '_ ` _ \
|  _| | |__| |___| | | |  __/ | | | | |
|_|   |_____\____|_| |_|\___|_| |_| |_|

Finite Element Method Solver for Chemical Engineering Applications
2D Finite Element Method Module - Source Code

Overview
--------
FEChem is a finite element method (FEM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.

Currently, only P1 elements (lin2, tri3, and quad4) and constant time steps are implemented.
Advective systems are stabilized via Streamline Upwind Petrov-Galerkin (SUPG).
The Navier-Stokes equation is stabilized via SUPG and Pressure-Stabilized Petrov-Galerkin (PSPG).

The source code is organized into the following folders:
- base      : core data structures (e.g., geometry, quadrature, scalars, vectors, writers)
- operator  : applies discretized weak form into global matrix equation
- shape     : shape functions and quadrature rules
- solver    : linear matrix solvers (e.g., LU, QR, GMRES)
- steady    : steady-state solvers
- transient : transient solvers

Base
----
- error    : defines errors returned by the library
- geom_bnd : stores 1D boundary geometry
- geom_dom : stores 2D domain geometry
- geom_itf : stores 1D interface geometry; interfaces join two domains
- itg_bnd  : stores quadrature data for a boundary
- itg_dom  : stores quadrature data for a domain
- itg_itf  : stores quadrature data for an interface
- mesh     : stores mesh nodes, elements, and physical regions
- read_gmsh: imports a mesh from a Gmsh file
- scl_bnd  : stores a scalar field defined on a boundary
- scl_dom  : stores a scalar field defined on a domain
- scl_itf  : stores a scalar field defined on an interface
- vars     : owns all geometry, quadrature, scalar, and vector data
- vec_bnd  : stores a 2D vector field defined on a boundary
- vec_dom  : stores a 2D vector field defined on a domain
- vec_itf  : stores a 2D vector field defined on an interface
- write_csv: writes scalar and vector values to CSV files
- write_vtu: writes scalar and vector values to VTU files

Operator
--------
Operator names use the form op + field type + location + term. For example,
opscl_dom_diff is a scalar domain diffusion operator. Files ending in _unity
are variants whose weighting coefficient is fixed at one.

- oper_base                    : defines the common operator trait and global matrix assembly helpers
- opscl_bnd_dir                : enforces a scalar Dirichlet boundary condition
- opscl_bnd_div                : adds prescribed boundary flow to the continuity equation
- opscl_bnd_neu                : applies a prescribed scalar normal flux
- opscl_bnd_out                : applies the weighted scalar advective outflow term
- opscl_bnd_out_unity          : applies the scalar advective outflow term with unit weight
- opscl_bnd_trn                : applies a scalar transfer (Robin) boundary condition
- opscl_dom_adv                : adds weighted advection to a scalar transport equation
- opscl_dom_adv_unity          : adds scalar advection with unit weight
- opscl_dom_den_time           : adds the backward-Euler density derivative to continuity
- opscl_dom_diff               : adds diffusion to a scalar transport equation
- opscl_dom_div                : adds velocity divergence to the continuity equation
- opscl_dom_pspg_steady        : adds steady PSPG stabilization to the continuity equation
- opscl_dom_pspg_time          : adds transient PSPG stabilization to the continuity equation
- opscl_dom_src                : adds a source to a scalar transport equation
- opscl_dom_supg_steady        : adds steady SUPG stabilization to weighted scalar transport
- opscl_dom_supg_steady_unity  : adds steady SUPG stabilization with unit weight
- opscl_dom_supg_time          : adds the transient SUPG time-residual term for weighted transport
- opscl_dom_supg_time_unity    : adds the transient SUPG time-residual term with unit weight
- opscl_dom_time               : adds a weighted backward-Euler scalar time derivative
- opscl_dom_time_unity         : adds a backward-Euler scalar time derivative with unit weight
- opscl_itf_cont               : enforces scalar value and flux continuity with a Lagrange multiplier
- opscl_itf_trn                : applies scalar transfer between two domains
- opvec_bnd_dir                : enforces a vector Dirichlet boundary condition
- opvec_bnd_pres               : applies a prescribed pressure boundary term to momentum transport
- opvec_dom_adv                : adds advection to the momentum equation
- opvec_dom_diff               : adds the viscous stress term to the momentum equation
- opvec_dom_pres               : adds the pressure-gradient term to the momentum equation
- opvec_dom_src                : adds a body-force source to the momentum equation
- opvec_dom_supg_steady        : adds steady SUPG stabilization to the momentum equation
- opvec_dom_supg_time          : adds the transient SUPG time-residual term to momentum transport
- opvec_dom_time               : adds a backward-Euler vector time derivative
- opvec_itf_cont               : enforces vector value and flux continuity with a Lagrange multiplier
- prelude                      : re-exports operator types for convenient internal imports

Shape
-----
- shape_base  : defines the common 1D and 2D shape-function traits
- shp1d_lin2  : implements a two-node linear line element
- shp2d_tri3  : implements a three-node linear triangular element
- shp2d_quad4 : implements a four-node bilinear quadrilateral element
- prelude     : re-exports shape types for convenient internal imports

Solver
------
- solver_base  : defines the common linear solver trait
- solver_gmres : solves sparse systems with the iterative GMRES method
- solver_lu    : solves sparse systems with an LU factorization
- solver_qr    : solves sparse systems with a QR factorization

Steady
------
- steady_base         : implements the common steady nonlinear solution loop
- steady_flow         : solves steady momentum transfer
- steady_heat         : solves steady heat transfer
- steady_heatflow     : solves coupled steady heat and momentum transfer
- steady_heatmass     : solves coupled steady heat and mass transfer
- steady_heatmassflow : solves coupled steady heat, mass, and momentum transfer
- steady_mass         : solves steady mass transfer
- steady_massflow     : solves coupled steady mass and momentum transfer

Transient
---------
- transient_base         : implements the common time-stepping and nonlinear solution loops
- transient_flow         : solves transient momentum transfer
- transient_heat         : solves transient heat transfer
- transient_heatflow     : solves coupled transient heat and momentum transfer
- transient_heatmass     : solves coupled transient heat and mass transfer
- transient_heatmassflow : solves coupled transient heat, mass, and momentum transfer
- transient_mass         : solves transient mass transfer
- transient_massflow     : solves coupled transient mass and momentum transfer
