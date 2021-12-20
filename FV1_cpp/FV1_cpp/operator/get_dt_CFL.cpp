#include "get_dt_CFL.h"

void get_dt_CFL
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	AssembledSolution&    assem_sol, 
	real&                 dx, 
	real&                 dt, 
	real&                 total_mass
)
{
	for (int i = 1; i < sim_params.cells + 1; i++)
	{
		if (assem_sol.h_BC[i] > solver_params.tol_dry)
		{
			real u = assem_sol.q_BC[i] / assem_sol.h_BC[i];
			real dtCFL = solver_params.CFL * dx / (abs(u) + sqrt(solver_params.g * assem_sol.h_BC[i]));
			dt = std::min(dt, dtCFL);
		}

		total_mass += assem_sol.h_BC[i] * dx;
	}
}