#include "malloc_fluxes.h"

void malloc_fluxes
(
	Fluxes&               fluxes,
	SimulationParameters& sim_params
)
{
	fluxes.mass     = new real[sim_params.cells + 1];
	fluxes.momentum = new real[sim_params.cells + 1];
}