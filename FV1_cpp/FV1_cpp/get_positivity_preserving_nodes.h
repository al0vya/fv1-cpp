#pragma once

#include <algorithm>

#include "SimulationParameters.h"
#include "AssembledSolution.h"
#include "SolverParameters.h"
#include "FaceValues.h"
#include "StarValues.h"

void get_positivity_preserving_nodes
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	FaceValues&           face_vals, 
	StarValues&           star_vals, 
	real*&                delta_west, 
	real*&                delta_east
);