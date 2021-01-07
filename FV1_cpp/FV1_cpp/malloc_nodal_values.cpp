#include "malloc_nodal_values.h"

void malloc_nodal_values
(
	NodalValues&          nodal_values,
	SimulationParameters& sim_params
)
{
	nodal_values.q = new real[sim_params.cells + 1];
	nodal_values.h = new real[sim_params.cells + 1];
	nodal_values.z = new real[sim_params.cells + 1];
	nodal_values.x = new real[sim_params.cells + 1];
}