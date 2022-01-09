#pragma once

#include <fstream>
#include <iostream>
#include <string>

#include "../classes/SimulationParameters.h"
#include "../classes/AssembledSolution.h"
#include "../classes/NodalValues.h"
#include "../classes/SaveInterval.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	NodalValues&          nodal_vals, 
	AssembledSolution&    assem_sol,
	SaveInterval&         saveint
);