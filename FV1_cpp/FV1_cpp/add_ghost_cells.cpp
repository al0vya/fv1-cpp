#include "add_ghost_cells.h"

void add_ghost_cells
(
	AssembledSolution&    assem_sol, 
	BoundaryConditions&   bcs, 
	SimulationParameters& sim_params)
{
	assem_sol.q_BC[0] = (bcs.q_imposed_up > 0) ? bcs.q_imposed_up : bcs.reflect_up * assem_sol.q_BC[1];
	assem_sol.q_BC[sim_params.cells + 1] = (bcs.q_imposed_down > 0) ? bcs.q_imposed_down: bcs.reflect_down* assem_sol.q_BC[sim_params.cells];

	assem_sol.h_BC[0] = (bcs.h_imposed_up > 0) ? bcs.h_imposed_up : assem_sol.h_BC[1];
	assem_sol.h_BC[sim_params.cells + 1] = (bcs.h_imposed_down > 0) ? bcs.h_imposed_down : assem_sol.h_BC[sim_params.cells];

	assem_sol.z_BC[0] = assem_sol.z_BC[1];
	assem_sol.z_BC[sim_params.cells + 1] = assem_sol.z_BC[sim_params.cells];
}