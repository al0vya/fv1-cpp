#include "get_nodal_values.h"

void get_nodal_values
(
	NodalValues&          nodal_vals,
	SimulationParameters& sim_params,
	BoundaryConditions&   bcs,
	real&                 dx,
	int                   test_case
)
{
	for (int i = 0; i < sim_params.cells + 1; i++) nodal_vals.x[i] = sim_params.xmin + i * dx;
	
	switch (test_case)
	{
	case 1:
	case 2:
	case 3:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_vals.z[i] = 0;
			nodal_vals.h[i] = h_init_overtop(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = (nodal_vals.x[i]) <= 32.5 ? bcs.ql : bcs.qr;
		}
		break;
	case 4:
	case 5:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_vals.z[i] = bed_data_c_property(nodal_vals.x[i]);
			nodal_vals.h[i] = h_init_c_property(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = (nodal_vals.x[i]) <= 32.5 ? bcs.ql : bcs.qr;
		}
		break;
	case 6:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_vals.z[i] = bed_data_c_property(nodal_vals.x[i]);
			nodal_vals.h[i] = h_init_overtop(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = (nodal_vals.x[i]) <= 32.5 ? bcs.ql : bcs.qr;
		}
		break;
	default:
		break;
	}
}