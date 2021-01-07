#include "malloc_assembled_solution.h"

void malloc_assembled_solution
(
	AssembledSolution&    assem_sol,
	SimulationParameters& sim_params
)
{
	assem_sol.qWithBC = new real[sim_params.cells + 2];
	assem_sol.hWithBC = new real[sim_params.cells + 2];
	assem_sol.zWithBC = new real[sim_params.cells + 2];
}