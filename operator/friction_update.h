#pragma once

#include "AssembledSolution.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "friction_implicit.h"

void friction_update
(
	AssembledSolution&    assem_sol, 
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	real&                 dt);