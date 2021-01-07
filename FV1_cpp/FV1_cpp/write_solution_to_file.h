#pragma once

#include <fstream>

#include "SimulationParameters.h"
#include "AssembledSolution.h"
#include "NodalValues.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	NodalValues&          nodal_vals, 
	AssembledSolution&    assem_sol
);