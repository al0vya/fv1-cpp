#pragma once

#include "../classes/AssembledSolution.h"
#include "../classes/BoundaryConditions.h"
#include "../classes/SimulationParameters.h"

void add_ghost_cells
(
	AssembledSolution&    assem_sol, 
	BoundaryConditions&   bcs, 
	SimulationParameters& sim_params);