# FV1-CPU

This is a model or *solver*, termed 'FV1-CPU', for the one-dimensional (1D) shallow water equations (SWE), which are a set of hyperbolic partial differential equations. FV1-CPU has three other counterpart solvers for the 1D SWE:

* [HFV1-CPU](github.com/al0vya/HFV1_cpp): an adaptive solver
* [FV1-GPU](github.com/al0vya/FV1_GPU): a parallelised solver
* HFV1-GPU: a parallelised adaptive solver

FV1-CPU is a finite volume (FV1) model. These models represent a physical domain of interest using a mesh, which comprises discrete elements or cells. A discretised form of the PDE is solved to obtain quantities of interest e.g. water height or velocity, over each cell. The discretised form can be obtained by simplifying differential or integral operators in the PDE into algebraic relations that can be more easily computed to solve the PDE. By solving over each and every cell, quantities of interest are obtained the over entire mesh, thereby modelling the physical domain.

