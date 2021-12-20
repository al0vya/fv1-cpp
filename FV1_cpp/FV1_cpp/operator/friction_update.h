#pragma once

#include "../classes/AssembledSolution.h"
#include "../classes/SimulationParameters.h"
#include "../classes/SolverParameters.h"
#include "friction_implicit.h"

void friction_update
(
	AssembledSolution&    assem_sol, 
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	real&                 dt);