#include "hInitialOvertopping.h"

real hInitialOvertopping(BoundaryConditions bcs, real z_int, real x_int)
{
	real etaWest = bcs.hl;
	real etaEast = bcs.hr;

	real h;

	h = (x_int <= 25) ? etaWest - z_int : (etaEast - z_int < 0) ? bcs.hr : etaEast - z_int;

	return h;
}