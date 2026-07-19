 _____ _____ ____ _                    
|  ___| ____/ ___| |__   ___ _ __ ___  
| |_  |  _|| |   | '_ \ / _ \ '_ ` _ \ 
|  _| | |__| |___| | | |  __/ | | | | |
|_|   |_____\____|_| |_|\___|_| |_| |_|

Finite Element Method Solver for Chemical Engineering Applications
2D Finite Element Method Module - Example Files

How to Use Examples
-------------------
First-time users are encouraged to run quickstart.rs

Change the directory to fechem_fem2d then run quickstart.rs:
`cargo run --release --example quickstart`

For the other examples, MOVE THE FILE to the "examples" before running:
`cargo run --release --example name_of_file`

You must MOVE THE FILE into the "examples" folder. Otherwise it will not run.
         ^^^^^^^^^^^^^
Correct : femchem_fem2d/examples/name_of_file.rs
WRONG   : femchem_fem2d/examples/physics/name_of_file.rs

List of Examples
----------------
Below is a list of examples and the concepts they introduce:

Heat Transfer
- quickstart.rs   | introduction to FEChem code
- heat_func.rs    | non-constant thermal properties
- heat_gmsh.rs    | importing custom meshes
- heat_multi.rs   | multiple domains
- heat_time.rs    | transient heat transfer simulation

Mass Transfer
- mass_react.rs   | steady-state diffusion-reaction system

Momentum Transfer
- flow_channel.rs | flow through a Z-shaped channel
- flow_lid.rs     | lid-driven cavity flow
- flow_vortex.rs  | vortex shedding past a cylinder

Heat and Momentum Transfer
- heatflow_conv.rs    | natural convection
- heatflow_heater.rs  | heat transfer from a solid to a moving fluid
- heatflow_lid.rs     | lid-driven cavity flow with heat transfer

Mass and Momentum Transfer
- massflow_lid.rs     | lid-driven cavity flow with a reaction
- massflow_channel.rs | reactive flow through a Z-shaped channel

Heat, Mass, and Momentum Transfer
- heatmassflow_lid.rs    | lid-driven cavity flow with a reaction and heat transfer
- heatmassflow_heater.rs | heat transfer from a solid to a moving and reacting fluid
