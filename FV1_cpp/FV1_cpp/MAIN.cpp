// Aliases
#include "real.h"

// File input/output
#include "write_solution_to_file.h"

// Solver steps
#include "get_nodal_values.h"
#include "get_modal_values.h"
#include "add_ghost_cells.h"
#include "friction_update.h"
#include "get_wet_dry_cells.h"
#include "get_face_values.h"
#include "get_positivity_preserving_nodes.h"
#include "fluxHLL.h"
#include "get_bar_values.h"
#include "fv1_operator.h"
#include "get_dt_CFL.h"

// Structures
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "BoundaryConditions.h"
#include "AssembledSolution.h"
#include "BarValues.h"
#include "Fluxes.h"
#include "FaceValues.h"
#include "NodalValues.h"
#include "StarValues.h"

// Solver/simulation setters
#include "set_simulation_parameters.h"
#include "set_solver_parameters.h"
#include "set_boundary_conditions.h"

int main
(
	int    argc,
	char** argv
)
{
	if (argc < 3)
	{
		printf
		(
			"Please run in the command line as follows:\n\n"
			"FV1_cpp.exe <TEST_CASE> <NUM_CELLS>"
		);

		exit(-1);
	}
	
	int test_case = strtol(argv[1], nullptr, 10); //set_test_case();
	int num_cells = strtol(argv[2], nullptr, 10); //set_num_cells();

	clock_t start = clock();
	
	// =========================================================== //
	// INITIALISATION OF VARIABLES AND INSTANTIATION OF STRUCTURES //
	// =========================================================== //

	// Structures
	SimulationParameters sim_params    = set_simulation_parameters(test_case, num_cells);
	SolverParameters     solver_params = set_solver_parameters();
	BoundaryConditions   bcs           = set_boundary_conditions(test_case);

	NodalValues       nodal_vals(num_cells + 1);
	AssembledSolution assem_sol (num_cells);
	FaceValues        face_vals (num_cells + 1);
	StarValues        star_vals (num_cells + 1);
	Fluxes            fluxes    (num_cells + 1);
	BarValues         bar_vals  (num_cells);

	int*  dry_cells  = new int [sim_params.cells + 2];
	real* eta_temp   = new real[sim_params.cells + 2];
	real* delta_west = new real[sim_params.cells + 1];
	real* delta_east = new real[sim_params.cells + 1];

	// Variables
	real dx       = (sim_params.xmax - sim_params.xmin) / sim_params.cells;
	real time_now = 0;
	real dt       = C(1e-3);

	// =========================================================== //
	
	get_nodal_values
	(
		nodal_vals, 
		sim_params, 
		bcs,
		dx,
		test_case
	);

	get_modal_values
	(
		assem_sol, 
		nodal_vals, 
		sim_params
	);

	while (time_now < sim_params.simulationTime)
	{
		time_now += dt;
		
		if (time_now - sim_params.simulationTime > 0)
		{
			time_now -= dt;
			dt = sim_params.simulationTime - time_now;
			time_now += dt;
		}

		add_ghost_cells
		(
			assem_sol, 
			bcs, 
			sim_params
		);

		if (sim_params.manning > 0) friction_update(assem_sol, sim_params, solver_params, dt);		

		get_wet_dry_cells
		(
			dry_cells, 
			sim_params, 
			assem_sol, 
			solver_params
		);

		get_face_values
		(
			sim_params, 
			assem_sol, 
			face_vals, 
			eta_temp
		);

		get_positivity_preserving_nodes
		(
			sim_params, 
			solver_params, 
			face_vals, 
			star_vals, 
			delta_west, 
			delta_east
		);

		fluxHLL
		(
			sim_params, 
			solver_params, 
			star_vals, 
			fluxes
		);

		get_bar_values
		(
			star_vals, 
			sim_params, 
			bar_vals
		);

		fv1_operator
		(
			sim_params, 
			dry_cells, 
			dx, 
			fluxes, 
			solver_params, 
			bar_vals, 
			assem_sol, 
			dt
		);

		// CFL time step adjustment 
		dt = 1e9;
		real total_mass = 0;

		get_dt_CFL
		(
			sim_params, 
			solver_params, 
			assem_sol, 
			dx, 
			dt, 
			total_mass
		);

		printf("Mass: %.17g, time step: %f, time: %f s\n", total_mass, dt, time_now);
	}

	write_solution_to_file
	(
		sim_params, 
		nodal_vals, 
		assem_sol
	);

	// =================== //
	// MEMORY DEALLOCATION //
	// =================== //

	delete[] dry_cells;
	delete[] eta_temp;
	delete[] delta_west;
	delete[] delta_east;

	// =================== //

	// print execution time to console 
	clock_t end = clock();

	real end_time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", end_time);

	return 0;
}