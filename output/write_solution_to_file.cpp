#include "write_solution_to_file.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	NodalValues&          nodal_vals, 
	AssembledSolution&    assem_sol,
	SaveInterval&         saveint
)
{
	std::string filename = "solution_data-" + std::to_string(saveint.count - 1) + ".csv";
	
	std::ofstream test;

	test.open(filename);

	test << "x,q,z,eta" << std::endl;

	for (int i = 0; i < sim_params.cells; i++)
	{
		test << (nodal_vals.x[i] + nodal_vals.x[i + 1]) / 2 << "," 
			 << assem_sol.q_BC[i + 1] << "," 
			 << assem_sol.z_BC[i + 1] << "," 
			 << std::fmax(assem_sol.z_BC[i + 1], assem_sol.h_BC[i + 1] + assem_sol.z_BC[i + 1]) 
			 << "\n";
	}

	test.close();
}