#pragma once
#include "../classes/real.h"

typedef struct BoundaryConditions
{
	real hl;
	real hr;

	real ql;
	real qr;

	real reflect_up;
	real reflect_down;

	real h_imposed_up;
	real q_imposed_up;

	real h_imposed_down;
	real q_imposed_down;

} BoundaryConditions;

/*
example of a designated initialiser
requires /std:c++latest
the (slight) advantage is bcs will be initialised upon construction
there is never a time when bcs is (partially) uninitialised

BoundaryConditions bcs = {
	.hl = C(2.0),
	.hr = C(0.0),
	.ql = C(0.0),
	.qr = C(0.0),
	.reflectUp = C(0.0),
	.reflectDown = C(0.0),
	.hImposedUp = C(0.0),
	.qxImposedUp = C(0.0),
	.hImposedDown = C(0.0),
	.qxImposedDown = C(0.0)
};
*/