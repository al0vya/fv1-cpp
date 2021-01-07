#include "free_assembled_solution.h"

void free_assembled_solution(AssembledSolution& assem_sol)
{
	delete[] assem_sol.qWithBC;
	delete[] assem_sol.hWithBC;
	delete[] assem_sol.zWithBC;
}