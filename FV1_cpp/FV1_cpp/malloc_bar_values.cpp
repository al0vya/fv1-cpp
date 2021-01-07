#include "malloc_bar_values.h"

void malloc_bar_values
(
	BarValues&             bar_vals,
	SimulationParameters& sim_params
)
{
	bar_vals.h = new real[sim_params.cells];
	bar_vals.z = new real[sim_params.cells];
}