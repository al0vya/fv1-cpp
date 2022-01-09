#pragma once

#include "AssembledSolution.h"
#include "BoundaryConditions.h"
#include "SimulationParameters.h"

void add_ghost_cells
(
	AssembledSolution&    assem_sol, 
	BoundaryConditions&   bcs, 
	SimulationParameters& sim_params);