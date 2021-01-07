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
#include "set_simulation_parameters_for_test_case.h"
#include "set_solver_parameters_for_test_case.h"
#include "set_boundary_conditions_for_test_case.h"
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

int main()
{
	int steps = 0;

	int test_case_selection = set_test_case();
	int num_cells = set_num_cells();

	clock_t start = clock();
	
	// =========================================================== //
	// INITIALISATION OF VARIABLES AND INSTANTIATION OF STRUCTURES //
	// =========================================================== //

	SimulationParameters sim_params    = set_simulation_parameters_for_test_case(test_case_selection, num_cells);
	SolverParameters     solver_params = set_solver_parameters_for_test_case();
	BoundaryConditions   bcs           = set_boundary_conditions_for_test_case(test_case_selection);

	NodalValues       nodal_values;
	AssembledSolution assem_sol;
	FaceValues        face_values;
	StarValues        star_values;
	Fluxes            fluxes;
	BarValues         bar_values;

	// =========================================================== //
	
	
	// coarsest cell size
	real dx = (sim_params.xmax - sim_params.xmin) / sim_params.cells;

	// allocate buffer for interfaces
	real* x_int = new real[sim_params.cells + 1];

	// initialise baseline mesh
	for (int i = 0; i < sim_params.cells + 1; i++) x_int[i] = sim_params.xmin + i * dx;

	// allocate buffers for flow nodes
	nodal_values.q = new real[sim_params.cells + 1];
	nodal_values.h = new real[sim_params.cells + 1];
	nodal_values.z = new real[sim_params.cells + 1];

	// initial interface values
	switch (test_case_selection)
	{
	case 1:
	case 2:
	case 3:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_values.z[i] = 0;
			nodal_values.h[i] = hInitialOvertopping(bcs, nodal_values.z[i], x_int[i]);
			nodal_values.q[i] = qInitial(bcs, x_int[i]);
		}
		break;
	case 4:
	case 5:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_values.z[i] = bedDataConservative(x_int[i]);
			nodal_values.h[i] = hInitialConservative(bcs, nodal_values.z[i], x_int[i]);
			nodal_values.q[i] = qInitial(bcs, x_int[i]);
		}
		break;
	case 6:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			nodal_values.z[i] = bedDataConservative(x_int[i]);
			nodal_values.h[i] = hInitialOvertopping(bcs, nodal_values.z[i], x_int[i]);
			nodal_values.q[i] = qInitial(bcs, x_int[i]);
		}
		break;
	default:
		break;
	}

	// allocate buffers for flow modes with ghost BCs
	assem_sol.qWithBC = new real[sim_params.cells + 2];
	assem_sol.hWithBC = new real[sim_params.cells + 2];
	assem_sol.zWithBC = new real[sim_params.cells + 2];

	// project to find the modes
	for (int i = 1; i < sim_params.cells + 1; i++)
	{
		assem_sol.qWithBC[i] = (nodal_values.q[i - 1] + nodal_values.q[i]) / 2;
		assem_sol.hWithBC[i] = (nodal_values.h[i - 1] + nodal_values.h[i]) / 2;
		assem_sol.zWithBC[i] = (nodal_values.z[i - 1] + nodal_values.z[i]) / 2;
	}

	// allocate true/false buffer for dry cells
	int* dry = new int[sim_params.cells + 2];

	real* etaTemp = new real[sim_params.cells + 2];

	// allocating buffers for eastern and western interface values
	face_values.q_east = new real[sim_params.cells + 1];
	face_values.h_east = new real[sim_params.cells + 1];
	face_values.eta_east = new real[sim_params.cells + 1];

	face_values.q_west = new real[sim_params.cells + 1];
	face_values.h_west = new real[sim_params.cells + 1];
	face_values.eta_west = new real[sim_params.cells + 1];

	// allocating buffers for positivity preserving nodes
	star_values.q_east = new real[sim_params.cells + 1];
	star_values.h_east = new real[sim_params.cells + 1];

	star_values.q_west = new real[sim_params.cells + 1];
	star_values.h_west = new real[sim_params.cells + 1];

	star_values.z_east = new real[sim_params.cells + 1];
	star_values.z_west = new real[sim_params.cells + 1];
	
	real* deltaWest = new real[sim_params.cells + 1];
	real* deltaEast = new real[sim_params.cells + 1];

	// allocating buffers for numerical fluxes from HLL solver
	real* massFlux = new real[sim_params.cells + 1];
	real* momentumFlux = new real[sim_params.cells + 1];

	// allocating buffers for positivity preserving MODES
	real* hBar = new real[sim_params.cells];
	real* zBar = new real[sim_params.cells];

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
			face_values.q_east[i] = assem_sol.qWithBC[i + 1];
			face_values.h_east[i] = assem_sol.hWithBC[i + 1];
			face_values.eta_east[i] = etaTemp[i + 1];

			face_values.q_west[i] = assem_sol.qWithBC[i];
			face_values.h_west[i] = assem_sol.hWithBC[i];
			face_values.eta_west[i] = etaTemp[i];
		}

		// correcting downwind and upwind eta values
		face_values.eta_east[sim_params.cells] = etaTemp[sim_params.cells] - assem_sol.hWithBC[sim_params.cells] + assem_sol.hWithBC[sim_params.cells + 1];
		face_values.eta_west[0] = etaTemp[1] - assem_sol.hWithBC[1] + assem_sol.hWithBC[0];

		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			real uWest = (face_values.h_west[i] <= solver_params.tol_dry) ? 0 : face_values.q_west[i] / face_values.h_west[i];
			real uEast = (face_values.h_east[i] <= solver_params.tol_dry) ? 0 : face_values.q_east[i] / face_values.h_east[i];

			real zWest = face_values.eta_west[i] - face_values.h_west[i];
			real zEast = face_values.eta_east[i] - face_values.h_east[i];

			real z_star_intermediate = max(zWest, zEast);

			deltaWest[i] = max(C(0.0), -(face_values.eta_west[i] - z_star_intermediate));
			deltaEast[i] = max(C(0.0), -(face_values.eta_east[i] - z_star_intermediate));

			star_values.h_west[i] = max(C(0.0), face_values.eta_west[i] - z_star_intermediate);
			star_values.q_west[i] = uWest * star_values.h_west[i];

			star_values.h_east[i] = max(C(0.0), face_values.eta_east[i] - z_star_intermediate);
			star_values.q_east[i] = uEast * star_values.h_east[i];

			star_values.z_west[i] = z_star_intermediate - deltaWest[i];
			star_values.z_east[i] = z_star_intermediate - deltaEast[i];
		}

		// initialising numerical fluxes
		fluxHLL(sim_params, solver_params, star_values.h_west, star_values.h_east, star_values.q_west, star_values.q_east, massFlux, momentumFlux);

		for (int i = 0; i < sim_params.cells; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			hBar[i] = (star_values.h_west[i + 1] + star_values.h_east[i]) / 2;

			// 1st order projection
			zBar[i] = (star_values.z_west[i + 1] - star_values.z_east[i]) / (2 * sqrt(C(3.0)));
		}

		// FV1 operator increment, skip ghosts cells
		for (int i = 1; i < sim_params.cells + 1; i++)
		{
			// skip increment in dry cells
			if (!dry[i])
			{
				real massIncrement = -(1 / dx) * (massFlux[i] - massFlux[i - 1]);
				real momentumIncrement = -(1 / dx) * (momentumFlux[i] - momentumFlux[i - 1] + 2 * sqrt(C(3.0)) * solver_params.g * hBar[i - 1] * zBar[i - 1]);

				real a = momentumFlux[i] - momentumFlux[i - 1];
				real b = 2 * sqrt(C(3.0)) * solver_params.g * hBar[i - 1] * zBar[i - 1];

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
		test << (x_int[i] + x_int[i + 1]) / 2 << "," << assem_sol.qWithBC[i + 1] << "," << assem_sol.zWithBC[i + 1] << "," << max(assem_sol.zWithBC[i + 1], assem_sol.hWithBC[i + 1] + assem_sol.zWithBC[i + 1]) << endl;
	}

	test.close();

	// delete buffers
	delete[] x_int;

	delete[] nodal_values.q;
	delete[] nodal_values.h;
	delete[] nodal_values.z;

	delete[] assem_sol.qWithBC;
	delete[] assem_sol.hWithBC;
	delete[] assem_sol.zWithBC;

	delete[] dry;

	delete[] etaTemp;

	delete[] face_values.q_east;
	delete[] face_values.h_east;
	delete[] face_values.eta_east;

	delete[] face_values.q_west;
	delete[] face_values.h_west;
	delete[] face_values.eta_west;

	delete[] star_values.q_west;
	delete[] star_values.h_west;

	delete[] star_values.q_east;
	delete[] star_values.h_east;

	delete[] star_values.z_west;
	delete[] star_values.z_east;

	delete[] deltaWest;
	delete[] deltaEast;

	delete[] massFlux;
	delete[] momentumFlux;

	delete[] hBar;
	delete[] zBar;

	// print execution time to console 
	clock_t end = clock();

	real end_time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", time);

	return 0;
}