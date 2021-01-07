#include "get_face_values.h"

void get_face_values
(
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol, 
	FaceValues&           face_vals, 
	real*&                eta_temp
)
{
	for (int i = 0; i < sim_params.cells + 2; i++) eta_temp[i] = assem_sol.h_BC[i] + assem_sol.z_BC[i];
	
	for (int i = 0; i < sim_params.cells + 1; i++)
	{
		face_vals.q_east[i] = assem_sol.q_BC[i + 1];
		face_vals.h_east[i] = assem_sol.h_BC[i + 1];
		face_vals.eta_east[i] = eta_temp[i + 1];

		face_vals.q_west[i] = assem_sol.q_BC[i];
		face_vals.h_west[i] = assem_sol.h_BC[i];
		face_vals.eta_west[i] = eta_temp[i];
	}

	// correcting downwind and upwind eta values
	face_vals.eta_east[sim_params.cells] = eta_temp[sim_params.cells] - assem_sol.h_BC[sim_params.cells] + assem_sol.h_BC[sim_params.cells + 1];
	face_vals.eta_west[0] = eta_temp[1] - assem_sol.h_BC[1] + assem_sol.h_BC[0];
}