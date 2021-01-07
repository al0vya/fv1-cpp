#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>

#include "real.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "BoundaryConditions.h"
#include "set_simulation_parameters.h"
#include "set_solver_parameters.h"
#include "set_boundary_conditions.h"
#include "set_num_cells.h"
#include "set_test_case.h"
#include "bedDataConservative.h"
#include "hInitialConservative.h"
#include "hInitialOvertopping.h"
#include "qInitial.h"
#include "frictionImplicit.h"
#include "uFaceValues.h"
#include "fluxHLL.h"
#include "AssembledSolution.h"
#include "BarValues.h"
#include "Fluxes.h"
#include "FaceValues.h"
#include "NodalValues.h"
#include "StarValues.h"

// Memory (de)allocators
#include "malloc_assembled_solution.h"
#include "malloc_bar_values.h"
#include "malloc_face_values.h"
#include "malloc_fluxes.h"
#include "malloc_nodal_values.h"
#include "malloc_star_values.h"
#include "free_assembled_solution.h"
#include "free_bar_values.h"
#include "free_face_values.h"
#include "free_fluxes.h"
#include "free_nodal_values.h"
#include "free_star_values.h"

int main()
{
	int steps = 0;

	int test_case_selection = set_test_case();
	int num_cells = set_num_cells();

	clock_t start = clock();
	
	// =========================================================== //
	// INITIALISATION OF VARIABLES AND INSTANTIATION OF STRUCTURES //
	// =========================================================== //

	// Structures
	SimulationParameters sim_params    = set_simulation_parameters(test_case_selection, num_cells);
	SolverParameters     solver_params = set_solver_parameters();
	BoundaryConditions   bcs           = set_boundary_conditions(test_case_selection);

	NodalValues       nodal_vals;
	AssembledSolution assem_sol;
	FaceValues        face_vals;
	StarValues        star_vals;
	Fluxes            fluxes;
	BarValues         bar_vals;

	// Memory allocation
	malloc_nodal_values(nodal_vals, sim_params);
	malloc_assembled_solution(assem_sol, sim_params);
	malloc_face_values(face_vals, sim_params);
	malloc_star_values(star_vals, sim_params);
	malloc_fluxes(fluxes, sim_params);
	malloc_bar_values(bar_vals, sim_params);

	int* dry = new int[sim_params.cells + 2];

	real* etaTemp = new real[sim_params.cells + 2];

	real* deltaWest = new real[sim_params.cells + 1];
	real* deltaEast = new real[sim_params.cells + 1];

	// =========================================================== //
	
	// coarsest cell size
	real dx = (sim_params.xmax - sim_params.xmin) / sim_params.cells;

	// initialise baseline mesh
	for (int i = 0; i < sim_params.cells + 1; i++) nodal_vals.x[i] = sim_params.xmin + i * dx;

	// initial interface values
	switch (test_case_selection)
	{
	case 1:
	case 2:
	case 3:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_vals.z[i] = 0;
			nodal_vals.h[i] = hInitialOvertopping(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = qInitial(bcs, nodal_vals.x[i]);
		}
		break;
	case 4:
	case 5:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_vals.z[i] = bedDataConservative(nodal_vals.x[i]);
			nodal_vals.h[i] = hInitialConservative(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = qInitial(bcs, nodal_vals.x[i]);
		}
		break;
	case 6:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_vals.z[i] = bedDataConservative(nodal_vals.x[i]);
			nodal_vals.h[i] = hInitialOvertopping(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = qInitial(bcs, nodal_vals.x[i]);
		}
		break;
	default:
		break;
	}

	// project to find the modes
	for (int i = 1; i < sim_params.cells + 1; i++)
	{
		assem_sol.qWithBC[i] = (nodal_vals.q[i - 1] + nodal_vals.q[i]) / 2;
		assem_sol.hWithBC[i] = (nodal_vals.h[i - 1] + nodal_vals.h[i]) / 2;
		assem_sol.zWithBC[i] = (nodal_vals.z[i - 1] + nodal_vals.z[i]) / 2;
	}	

	real timeNow = 0;
	real dt = C(1e-3);

	int firstTimeStep = 1;

	// file stream for recording cumulative clock time vs sim time
	ofstream data;

	data.open("clock_time_vs_sim_time.csv");

	data.precision(15);

	data << "sim_time,clock_time" << std::endl;

	while (timeNow < sim_params.simulationTime)
	{
		timeNow += dt;
		
		if (timeNow - sim_params.simulationTime > 0)
		{
			timeNow -= dt;
			dt = sim_params.simulationTime - timeNow;
			timeNow += dt;
		}

		// adding ghost boundary conditions
		assem_sol.qWithBC[0] = (bcs.qxImposedUp > 0) ? bcs.qxImposedUp : bcs.reflectUp * assem_sol.qWithBC[1];
		assem_sol.qWithBC[sim_params.cells + 1] = (bcs.qxImposedDown > 0) ? bcs.qxImposedDown : bcs.reflectDown * assem_sol.qWithBC[sim_params.cells];

		assem_sol.hWithBC[0] = (bcs.hImposedUp > 0) ? bcs.hImposedUp : assem_sol.hWithBC[1];
		assem_sol.hWithBC[sim_params.cells + 1] = (bcs.hImposedDown > 0) ? bcs.hImposedDown : assem_sol.hWithBC[sim_params.cells];

		assem_sol.zWithBC[0] = assem_sol.zWithBC[1];
		assem_sol.zWithBC[sim_params.cells + 1] = assem_sol.zWithBC[sim_params.cells];

		if (sim_params.manning > 0)
		{
			for (int i = 1; i < sim_params.cells + 1; i++)
			{
				assem_sol.qWithBC[i] += frictionImplicit(sim_params, solver_params, dt, assem_sol.hWithBC[i], assem_sol.qWithBC[i]);
			}
		}

		// ghost cells are always 0 i.e. false/wet
		dry[0] = false;
		dry[sim_params.cells + 1] = false;

		// initialising dry vs wet cells, ignore ghost cells
		for (int i = 1; i < sim_params.cells + 1; i++)
		{
			real hLocal = assem_sol.hWithBC[i];
			real hBackward = assem_sol.hWithBC[i - 1];
			real hForward = assem_sol.hWithBC[i + 1];

			real hMax = max(hLocal, hBackward);
			hMax = max(hForward, hMax);

			dry[i] = (hMax <= solver_params.tol_dry);
		}

		for (int i = 0; i < sim_params.cells + 2; i++)
		{
			etaTemp[i] = assem_sol.hWithBC[i] + assem_sol.zWithBC[i];
		}

		// initialising interface values
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			face_vals.q_east[i] = assem_sol.qWithBC[i + 1];
			face_vals.h_east[i] = assem_sol.hWithBC[i + 1];
			face_vals.eta_east[i] = etaTemp[i + 1];

			face_vals.q_west[i] = assem_sol.qWithBC[i];
			face_vals.h_west[i] = assem_sol.hWithBC[i];
			face_vals.eta_west[i] = etaTemp[i];
		}

		// correcting downwind and upwind eta values
		face_vals.eta_east[sim_params.cells] = etaTemp[sim_params.cells] - assem_sol.hWithBC[sim_params.cells] + assem_sol.hWithBC[sim_params.cells + 1];
		face_vals.eta_west[0] = etaTemp[1] - assem_sol.hWithBC[1] + assem_sol.hWithBC[0];

		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			real uWest = (face_vals.h_west[i] <= solver_params.tol_dry) ? 0 : face_vals.q_west[i] / face_vals.h_west[i];
			real uEast = (face_vals.h_east[i] <= solver_params.tol_dry) ? 0 : face_vals.q_east[i] / face_vals.h_east[i];

			real zWest = face_vals.eta_west[i] - face_vals.h_west[i];
			real zEast = face_vals.eta_east[i] - face_vals.h_east[i];

			real z_star_intermediate = max(zWest, zEast);

			deltaWest[i] = max(C(0.0), -(face_vals.eta_west[i] - z_star_intermediate));
			deltaEast[i] = max(C(0.0), -(face_vals.eta_east[i] - z_star_intermediate));

			star_vals.h_west[i] = max(C(0.0), face_vals.eta_west[i] - z_star_intermediate);
			star_vals.q_west[i] = uWest * star_vals.h_west[i];

			star_vals.h_east[i] = max(C(0.0), face_vals.eta_east[i] - z_star_intermediate);
			star_vals.q_east[i] = uEast * star_vals.h_east[i];

			star_vals.z_west[i] = z_star_intermediate - deltaWest[i];
			star_vals.z_east[i] = z_star_intermediate - deltaEast[i];
		}

		// initialising numerical fluxes
		fluxHLL(sim_params, solver_params, star_vals.h_west, star_vals.h_east, star_vals.q_west, star_vals.q_east, fluxes.mass, fluxes.momentum);

		for (int i = 0; i < sim_params.cells; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			bar_vals.h[i] = (star_vals.h_west[i + 1] + star_vals.h_east[i]) / 2;

			// 1st order projection
			bar_vals.z[i] = (star_vals.z_west[i + 1] - star_vals.z_east[i]) / (2 * sqrt(C(3.0)));
		}

		// FV1 operator increment, skip ghosts cells
		for (int i = 1; i < sim_params.cells + 1; i++)
		{
			// skip increment in dry cells
			if (!dry[i])
			{
				real massIncrement = -(1 / dx) * (fluxes.mass[i] - fluxes.mass[i - 1]);
				real momentumIncrement = -(1 / dx) * (fluxes.momentum[i] - fluxes.momentum[i - 1] + 2 * sqrt(C(3.0)) * solver_params.g * bar_vals.h[i - 1] * bar_vals.z[i - 1]);

				real a = fluxes.momentum[i] - fluxes.momentum[i - 1];
				real b = 2 * sqrt(C(3.0)) * solver_params.g * bar_vals.h[i - 1] * bar_vals.z[i - 1];

				assem_sol.hWithBC[i] += dt * massIncrement;
				assem_sol.qWithBC[i] = (assem_sol.hWithBC[i] <= solver_params.tol_dry) ? 0 : assem_sol.qWithBC[i] + dt * momentumIncrement;
			}
		}

		// CFL time step adjustmenttotal_mass 
		dt = 1e9;
		real total_mass = 0;

		for (int i = 1; i < sim_params.cells + 1; i++)
		{
			if (assem_sol.hWithBC[i] > solver_params.tol_dry)
			{
				real u = assem_sol.qWithBC[i] / assem_sol.hWithBC[i];
				real dtCFL = solver_params.CFL * dx / (abs(u) + sqrt(solver_params.g * assem_sol.hWithBC[i]));
				dt = min(dt, dtCFL);
			}

			total_mass  += assem_sol.hWithBC[i] * dx;
		}

		steps++;

		// recording cumulative clock time against sim time
		if (steps % 100 == 0)
		{
			clock_t current = clock();
			real time = (real)(current - start) / CLOCKS_PER_SEC * C(1000.0);

			data << timeNow << "," << time << endl;
		}

		printf("Mass: %.17g, time step: %f, time: %f s\n", total_mass , dt, timeNow);
	}

	
	// ensures recording of final clock time if steps % 100 != 0
	clock_t current = clock();
	real time = (real)(current - start) / CLOCKS_PER_SEC * C(1000.0);

	data << timeNow << "," << time << endl;

	data.close();

	// seperate file stream for x, q, max(eta, h) and z
	ofstream test;

	test.open("FV1_data.csv");

	test << "x,q,z,eta" << endl;

	for (int i = 0; i < sim_params.cells; i++)
	{
		test << (nodal_vals.x[i] + nodal_vals.x[i + 1]) / 2 << "," << assem_sol.qWithBC[i + 1] << "," << assem_sol.zWithBC[i + 1] << "," << max(assem_sol.zWithBC[i + 1], assem_sol.hWithBC[i + 1] + assem_sol.zWithBC[i + 1]) << endl;
	}

	test.close();

	// delete buffers
	free_nodal_values(nodal_vals);
	free_assembled_solution(assem_sol);
	free_face_values(face_vals);
	free_star_values(star_vals);
	free_bar_values(bar_vals);
	free_fluxes(fluxes);

	delete[] dry;
	delete[] etaTemp;
	delete[] deltaWest;
	delete[] deltaEast;

	// print execution time to console 
	clock_t end = clock();

	real end_time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", time);

	return 0;
}