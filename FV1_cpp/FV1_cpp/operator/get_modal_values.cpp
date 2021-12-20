#include "get_modal_values.h"

void get_modal_values
(
	AssembledSolution&    assem_sol, 
	NodalValues&          nodal_vals, 
	SimulationParameters& sim_params
)
{
	for (int i = 1; i < sim_params.cells + 1; i++)
	{
		assem_sol.q_BC[i] = (nodal_vals.q[i - 1] + nodal_vals.q[i]) / 2;
		assem_sol.h_BC[i] = (nodal_vals.h[i - 1] + nodal_vals.h[i]) / 2;
		assem_sol.z_BC[i] = (nodal_vals.z[i - 1] + nodal_vals.z[i]) / 2;
	}
}