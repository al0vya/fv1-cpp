#pragma once

#include "SimulationParameters.h"
#include "NodalValues.h"
#include "BoundaryConditions.h"
#include "h_init_c_property.h"
#include "h_init_overtop.h"
#include "bed_data_c_property.h"

void get_nodal_values
(
	NodalValues&          nodal_vals,
	SimulationParameters& sim_params,
	BoundaryConditions&   bcs,
	real&                 dx,
	int                   test_case
);