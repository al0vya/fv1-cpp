#pragma once

#include "../classes/AssembledSolution.h"
#include "../classes/SimulationParameters.h"
#include "../classes/FaceValues.h"

void get_face_values
(
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol, 
	FaceValues&           face_vals, 
	real*&                eta_temp
);