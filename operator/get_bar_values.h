#pragma once

#include <algorithm>

#include "../classes/StarValues.h"
#include "../classes/SimulationParameters.h"
#include "../classes/BarValues.h"

void get_bar_values
(
	StarValues&           star_vals, 
	SimulationParameters& sim_params, 
	BarValues&            bar_vals
);