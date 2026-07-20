 _____ _____ ____ _                    
|  ___| ____/ ___| |__   ___ _ __ ___  
| |_  |  _|| |   | '_ \ / _ \ '_ ` _ \ 
|  _| | |__| |___| | | |  __/ | | | | |
|_|   |_____\____|_| |_|\___|_| |_| |_|

Finite Element Method Solver for Chemical Engineering Applications
2D Finite Element Method Module
Copyright (c) 2026 FEChem Development Team

Overview
--------
FEChem is a finite element method (FEM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.
It is a pure Rust library. Simulations are configured as Rust code, and no separate installation or build system is required.

Note: The project is under active development. The public API and numerical models may change between releases.

Features
--------
- Steady or transient simulations of heat, mass, or momentum transport and couplings thereof
- Non-constant inputs, such as temperature and concentration-dependent properties and reactions
- Multiple domains, boundaries, and interfaces, which can be imported from Gmsh meshes
- Streamline upwind Petrov-Galerkin (SUPG) stabilization for advective transport
- SUPG and Pressure-Stabilized Petrov-Galerkin (PSPG) stabilization for the Navier-Stokes equation
- Direct LU and QR linear solvers and an iterative GMRES solver
- Outputs simulation results as CSV or VTU files for post-processing

Requirements
------------
- A Rust toolchain with Rust 2024 edition support
- Gmsh for custom meshes (optional)
- ParaView for visualizing output (optional)

Quick Start
-----------
Clone or download the project, change to the fechem_fem2d directory, and run `cargo run --release --example quickstart`
This example solves a steady heat transfer problem on a square mesh and writes the temperature field to `examples/output_quickstart/temp_0.vtu`.
Other examples can be found in the `examples` folder. An overview of the examples is provided in `examples/readme.txt`.

Project Structure
-----------------
A guide to the source code is provided in `src/readme.txt`.

Contributing
------------
Due to limited maintenance capacity, we are not currently accepting external contributions.
However, users are welcome to fork and modify the repository in accordance with the MIT License.

License
-------
FEChem is available under the MIT License. See license.txt for details.
