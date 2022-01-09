#pragma once

#include <cmath>

#include "SimulationParameters.h"
#include "AssembledSolution.h"
#include "SolverParameters.h"

void get_wet_dry_cells
(
	int*                  dry, 
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol, 
	SolverParameters&     solver_params
);