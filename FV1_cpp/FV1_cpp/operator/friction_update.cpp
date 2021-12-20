#include "friction_update.h"

void friction_update
(
	AssembledSolution&    assem_sol, 
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	real&                 dt)
{
	for (int i = 1; i < sim_params.cells + 1; i++)
	{
		assem_sol.q_BC[i] += friction_implicit(sim_params, solver_params, dt, assem_sol.h_BC[i], assem_sol.q_BC[i]);
	}
}