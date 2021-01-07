#include "free_fluxes.h"

void free_fluxes(Fluxes& fluxes)
{
	delete[] fluxes.mass;
	delete[] fluxes.momentum;
}