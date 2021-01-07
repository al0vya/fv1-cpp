#include "malloc_assembled_solution.h"

void malloc_assembled_solution
(
	AssembledSolution&    assem_sol,
	SimulationParameters& sim_params
)
{
	assem_sol.q_BC = new real[sim_params.cells + 2];
	assem_sol.h_BC = new real[sim_params.cells + 2];
	assem_sol.z_BC= new real[sim_params.cells + 2];
}