# FV1-CPU

## Model description

This is a shallow water model or *solver*, termed 'FV1-CPU', for the one-dimensional (1D) [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations) (SWE), which are a set of hyperbolic partial differential equations (PDEs). FV1-CPU has three other counterpart solvers for the 1D SWE:

* [HFV1-CPU](github.com/al0vya/HFV1_cpp): an adaptive solver
* [FV1-GPU](github.com/al0vya/FV1_GPU): a parallelised solver
* HFV1-GPU: a parallelised adaptive solver

FV1-CPU is a finite volume (FV1) model. These models represent a physical domain of interest using a mesh, which comprises discrete elements or cells. A discretised form of the PDE is solved to obtain quantities of interest e.g. water height or velocity, over each cell. The discretised form can be obtained by simplifying differential or integral operators in the PDE into algebraic relations that can be more easily computed to solve the PDE. By solving over each and every cell, quantities of interest are obtained over the entire mesh, thereby modelling the physical domain. An example of a physical domain that can be modelled using the 1D SWE is a channel (colloquially, a river).

## Running the model

Six test case simulations can be run using FV1-CPU to simulate different situations that can arise during shallow water flow in a channel, shown below. After building and then running the executable, the user must select which test case to run and how many cells are to comprise the mesh.

To understand the following test case animations, imagine viewing the channel from the side on i.e. a cross sectional view.

### Wet dam break

<img src="https://github.com/al0vya/FV1_GPU/blob/master/FV1_GPU_1D/test_case_gifs/wet_dam_break.gif" width="50%" height="50%">

### Dry dam break

<img src="https://github.com/al0vya/FV1_GPU/blob/master/FV1_GPU_1D/test_case_gifs/dry_dam_break.gif" width="50%" height="50%">

### Dry dam break with friction

<img src="https://github.com/al0vya/FV1_GPU/blob/master/FV1_GPU_1D/test_case_gifs/dry_dam_break_fric.gif" width="50%" height="50%">

### Water at rest

<img src="https://github.com/al0vya/FV1_GPU/blob/master/FV1_GPU_1D/test_case_gifs/wet_c_property.gif" width="50%" height="50%">

### Water at rest with dry zones

<img src="https://github.com/al0vya/FV1_GPU/blob/master/FV1_GPU_1D/test_case_gifs/wet_dry_c_property.gif" width="50%" height="50%">

### Dam break with water overtopping a building

<img src="https://github.com/al0vya/FV1_GPU/blob/master/FV1_GPU_1D/test_case_gifs/building_overtopping.gif" width="50%" height="50%">
