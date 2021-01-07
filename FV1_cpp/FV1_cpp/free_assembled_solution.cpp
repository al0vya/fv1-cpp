#include "free_assembled_solution.h"

void free_assembled_solution(AssembledSolution& assem_sol)
{
	delete[] assem_sol.q_BC;
	delete[] assem_sol.h_BC;
	delete[] assem_sol.z_BC;
}