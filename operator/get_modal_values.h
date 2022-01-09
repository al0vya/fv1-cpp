#pragma once

#include "AssembledSolution.h"
#include "NodalValues.h"
#include "SimulationParameters.h"

void get_modal_values
(
	AssembledSolution&    assem_sol, 
	NodalValues&          nodal_vals, 
	SimulationParameters& sim_params
);