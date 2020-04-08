#include "qInitial.h"

real qInitial(BoundaryConditions bcs, real x_int)
{
	if (x_int <= 32.5)
	{
		return bcs.ql;
	}
	else
	{
		return bcs.qr;
	}
}