#pragma once

#include "../classes/AssembledSolution.h"
#include "../classes/NodalValues.h"
#include "../classes/SimulationParameters.h"

void get_modal_values
(
	AssembledSolution&    assem_sol, 
	NodalValues&          nodal_vals, 
	SimulationParameters& sim_params
);