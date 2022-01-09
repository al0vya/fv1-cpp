#pragma once

#include <algorithm>

#include "../classes/SimulationParameters.h"
#include "../classes/SolverParameters.h"
#include "../classes/Fluxes.h"
#include "../classes/AssembledSolution.h"
#include "../classes/BarValues.h"

void fv1_operator
(
	SimulationParameters& sim_params, 
	int*&                 dry, 
	real&                 dx, 
	Fluxes&               fluxes, 
	SolverParameters&     solver_params, 
	BarValues&            bar_vals, 
	AssembledSolution&    assem_sol, 
	real&                 dt
);