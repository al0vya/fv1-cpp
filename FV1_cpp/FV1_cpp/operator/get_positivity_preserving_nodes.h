#pragma once

#include <algorithm>

#include "../classes/SimulationParameters.h"
#include "../classes/AssembledSolution.h"
#include "../classes/SolverParameters.h"
#include "../classes/FaceValues.h"
#include "../classes/StarValues.h"

void get_positivity_preserving_nodes
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	FaceValues&           face_vals, 
	StarValues&           star_vals, 
	real*&                delta_west, 
	real*&                delta_east
);