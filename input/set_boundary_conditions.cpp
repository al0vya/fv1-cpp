#include "set_boundary_conditions.h"
#include "real.h"

BoundaryConditions set_boundary_conditions(int test_case_selection)
{
	BoundaryConditions bcs;

	switch (test_case_selection)
	{
	case 1:
		bcs.hl = C(6.0);
		bcs.hr = C(2.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflect_up = C(1.0);
		bcs.reflect_down = C(1.0);
		bcs.h_imposed_up = C(0.0);
		bcs.q_imposed_up = C(0.0);
		bcs.h_imposed_down = C(0.0);
		bcs.q_imposed_down = C(0.0);
		break;
	case 2:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflect_up = C(1.0);
		bcs.reflect_down = C(1.0);
		bcs.h_imposed_up = C(0.0);
		bcs.q_imposed_up = C(0.0);
		bcs.h_imposed_down = C(0.0);
		bcs.q_imposed_down = C(0.0);
	case 3:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflect_up = C(1.0);
		bcs.reflect_down = C(1.0);
		bcs.h_imposed_up = C(0.0);
		bcs.q_imposed_up = C(0.0);
		bcs.h_imposed_down = C(0.0);
		bcs.q_imposed_down = C(0.0);
		break;
	case 4:
		bcs.hl = C(6.0);
		bcs.hr = C(6.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflect_up = C(1.0);
		bcs.reflect_down = C(1.0);
		bcs.h_imposed_up = C(0.0);
		bcs.q_imposed_up = C(0.0);
		bcs.h_imposed_down = C(0.0);
		bcs.q_imposed_down = C(0.0);
		break;
	case 5:
		bcs.hl = C(2.0);
		bcs.hr = C(2.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflect_up = C(1.0);
		bcs.reflect_down = C(1.0);
		bcs.h_imposed_up = C(0.0);
		bcs.q_imposed_up = C(0.0);
		bcs.h_imposed_down = C(0.0);
		bcs.q_imposed_down = C(0.0);
		break;
	case 6:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		bcs.ql = C(0.0);
		bcs.qr = C(0.0);
		bcs.reflect_up = C(1.0);
		bcs.reflect_down = C(1.0);
		bcs.h_imposed_up = C(0.0);
		bcs.q_imposed_up = C(0.0);
		bcs.h_imposed_down = C(0.0);
		bcs.q_imposed_down = C(0.0);
		break;
	default:
		break;
	}

	return bcs;
}