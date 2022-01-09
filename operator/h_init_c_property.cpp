#include "h_init_c_property.h"

real h_init_c_property(BoundaryConditions bcs, real z_int, real x_int)
{
	real etaWest = bcs.hl;
	real etaEast = bcs.hr;

	real h = etaWest - z_int;

	return (x_int <= 25) ? ((h < 0) ? bcs.hl : h) : etaEast - z_int;
}