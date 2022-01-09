#pragma once

#include <cmath>

#include "AssembledSolution.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"

void get_dt_CFL
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	AssembledSolution&    assem_sol, 
	real&                 dx, 
	real&                 dt, 
	real&                 total_mass
);