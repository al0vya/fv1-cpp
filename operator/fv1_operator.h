#pragma once

#include <cmath>

#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "Fluxes.h"
#include "AssembledSolution.h"
#include "BarValues.h"

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