#include "malloc_star_values.h"

void malloc_star_values
(
	StarValues&           star_values,
	SimulationParameters& sim_params
)
{
	star_values.q_east = new real[sim_params.cells + 1];
	star_values.h_east = new real[sim_params.cells + 1];

	star_values.q_west = new real[sim_params.cells + 1];
	star_values.h_west = new real[sim_params.cells + 1];

	star_values.z_east = new real[sim_params.cells + 1];
	star_values.z_west = new real[sim_params.cells + 1];
}