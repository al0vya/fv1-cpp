#include "hInitial.h"

real hInitial(BoundaryConditions bcs, real z_int, real x_int)
{
	real etaLeft = bcs.hl;
	real etaRight = bcs.hr;
	real a;

	if (x_int <= 25)
	{
		a = etaLeft - z_int;

		if (a < 0)
		{
			return 0;
		}
		else
		{
			return a;
		}
	}
	else
	{
		a = etaRight - z_int;

		if (a < 0)
		{
			return 0;
		}
		else
		{
			return a;
		}
	}
}