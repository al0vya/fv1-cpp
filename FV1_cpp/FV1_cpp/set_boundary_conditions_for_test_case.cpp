#include "set_boundary_conditions_for_test_case.h"
#include "real.h"

BoundaryConditions set_boundary_conditions_for_test_case(int test_case_selection)
{
	BoundaryConditions bcs;

	switch (test_case_selection)
	{
	case 1:
		bcs.hl = C(6.0);
		bcs.hr = C(2.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflectUp = C(1.0);
		bcs.reflectDown = C(1.0);
		bcs.hImposedUp = C(0.0);
		bcs.qxImposedUp = C(0.0);
		bcs.hImposedDown = C(0.0);
		bcs.qxImposedDown = C(0.0);
		break;
	case 2:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflectUp = C(1.0);
		bcs.reflectDown = C(1.0);
		bcs.hImposedUp = C(0.0);
		bcs.qxImposedUp = C(0.0);
		bcs.hImposedDown = C(0.0);
		bcs.qxImposedDown = C(0.0);
	case 3:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflectUp = C(1.0);
		bcs.reflectDown = C(1.0);
		bcs.hImposedUp = C(0.0);
		bcs.qxImposedUp = C(0.0);
		bcs.hImposedDown = C(0.0);
		bcs.qxImposedDown = C(0.0);
		break;
	case 4:
		bcs.hl = C(6.0);
		bcs.hr = C(6.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflectUp = C(1.0);
		bcs.reflectDown = C(1.0);
		bcs.hImposedUp = C(0.0);
		bcs.qxImposedUp = C(0.0);
		bcs.hImposedDown = C(0.0);
		bcs.qxImposedDown = C(0.0);
		break;
	case 5:
		bcs.hl = C(2.0);
		bcs.hr = C(2.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflectUp = C(1.0);
		bcs.reflectDown = C(1.0);
		bcs.hImposedUp = C(0.0);
		bcs.qxImposedUp = C(0.0);
		bcs.hImposedDown = C(0.0);
		bcs.qxImposedDown = C(0.0);
		break;
	case 6:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflectUp = C(1.0);
		bcs.reflectDown = C(1.0);
		bcs.hImposedUp = C(0.0);
		bcs.qxImposedUp = C(0.0);
		bcs.hImposedDown = C(0.0);
		bcs.qxImposedDown = C(0.0);
		break;
	default:
		break;
	}

	return bcs;
}