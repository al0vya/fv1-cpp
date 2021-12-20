#pragma once

#include <algorithm>

#include "../classes/SimulationParameters.h"
#include "../classes/AssembledSolution.h"
#include "../classes/SolverParameters.h"

void get_wet_dry_cells
(
	int*                  dry, 
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol, 
	SolverParameters&     solver_params
);