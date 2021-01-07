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

int main()
{
	clock_t start = clock();
	int steps = 0;

	int test_case_selection = set_test_case();
	int num_cells = set_num_cells();

	SimulationParameters sim_params = set_simulation_parameters_for_test_case(test_case_selection, num_cells);
	SolverParameters solver_params = set_solver_parameters_for_test_case();
	BoundaryConditions bcs = set_boundary_conditions_for_test_case(test_case_selection);

	// coarsest cell size
	real dx = (sim_params.xmax - sim_params.xmin) / sim_params.cells;

	// allocate buffer for interfaces
	real* xInt = new real[sim_params.cells + 1];

	// initialise baseline mesh
	for (int i = 0; i < sim_params.cells + 1; i++) xInt[i] = sim_params.xmin + i * dx;

	// allocate buffers for flow nodes
	real* qInt = new real[sim_params.cells + 1];
	real* hInt = new real[sim_params.cells + 1];
	real* zInt = new real[sim_params.cells + 1];

	// initial interface values
	switch (test_case_selection)
	{
	case 1:
	case 2:
	case 3:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			zInt[i] = 0;
			hInt[i] = hInitialOvertopping(bcs, zInt[i], xInt[i]);
			qInt[i] = qInitial(bcs, xInt[i]);
		}
		break;
	case 4:
	case 5:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			zInt[i] = bedDataConservative(xInt[i]);
			hInt[i] = hInitialConservative(bcs, zInt[i], xInt[i]);
			qInt[i] = qInitial(bcs, xInt[i]);
		}
		break;
	case 6:
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			zInt[i] = bedDataConservative(xInt[i]);
			hInt[i] = hInitialOvertopping(bcs, zInt[i], xInt[i]);
			qInt[i] = qInitial(bcs, xInt[i]);
		}
		break;
	default:
		break;
	}

	// allocate buffers for flow modes with ghost BCs
	real* qWithBC = new real[sim_params.cells + 2];
	real* hWithBC = new real[sim_params.cells + 2];
	real* zWithBC = new real[sim_params.cells + 2];

	// project to find the modes
	for (int i = 1; i < sim_params.cells + 1; i++)
	{
		qWithBC[i] = (qInt[i - 1] + qInt[i]) / 2;
		hWithBC[i] = (hInt[i - 1] + hInt[i]) / 2;
		zWithBC[i] = (zInt[i - 1] + zInt[i]) / 2;
	}

	// allocate true/false buffer for dry cells
	int* dry = new int[sim_params.cells + 2];

	real* etaTemp = new real[sim_params.cells + 2];

	// allocating buffers for eastern and western interface values
	real* qEast = new real[sim_params.cells + 1];
	real* hEast = new real[sim_params.cells + 1];
	real* etaEast = new real[sim_params.cells + 1];

	real* qWest = new real[sim_params.cells + 1];
	real* hWest = new real[sim_params.cells + 1];
	real* etaWest = new real[sim_params.cells + 1];

	// allocating buffers for positivity preserving nodes
	real* qEastStar = new real[sim_params.cells + 1];
	real* hEastStar = new real[sim_params.cells + 1];

	real* qWestStar = new real[sim_params.cells + 1];
	real* hWestStar = new real[sim_params.cells + 1];

	real* zWestStar = new real[sim_params.cells + 1];
	real* zEastStar = new real[sim_params.cells + 1];
	
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
		qWithBC[0] = (bcs.qxImposedUp > 0) ? bcs.qxImposedUp : bcs.reflectUp * qWithBC[1];
		qWithBC[sim_params.cells + 1] = (bcs.qxImposedDown > 0) ? bcs.qxImposedDown : bcs.reflectDown * qWithBC[sim_params.cells];

		hWithBC[0] = (bcs.hImposedUp > 0) ? bcs.hImposedUp : hWithBC[1];
		hWithBC[sim_params.cells + 1] = (bcs.hImposedDown > 0) ? bcs.hImposedDown : hWithBC[sim_params.cells];

		zWithBC[0] = zWithBC[1];
		zWithBC[sim_params.cells + 1] = zWithBC[sim_params.cells];

		if (sim_params.manning > 0)
		{
			for (int i = 1; i < sim_params.cells + 1; i++)
			{
				qWithBC[i] += frictionImplicit(sim_params, solver_params, dt, hWithBC[i], qWithBC[i]);
			}
		}

		// ghost cells are always 0 i.e. false/wet
		dry[0] = false;
		dry[sim_params.cells + 1] = false;

		// initialising dry vs wet cells, ignore ghost cells
		for (int i = 1; i < sim_params.cells + 1; i++)
		{
			real hLocal = hWithBC[i];
			real hBackward = hWithBC[i - 1];
			real hForward = hWithBC[i + 1];

			real hMax = max(hLocal, hBackward);
			hMax = max(hForward, hMax);

			dry[i] = (hMax <= solver_params.tol_dry);
		}

		for (int i = 0; i < sim_params.cells + 2; i++)
		{
			etaTemp[i] = hWithBC[i] + zWithBC[i];
		}

		// initialising interface values
		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			qEast[i] = qWithBC[i + 1];
			hEast[i] = hWithBC[i + 1];
			etaEast[i] = etaTemp[i + 1];

			qWest[i] = qWithBC[i];
			hWest[i] = hWithBC[i];
			etaWest[i] = etaTemp[i];
		}

		// correcting downwind and upwind eta values
		etaEast[sim_params.cells] = etaTemp[sim_params.cells] - hWithBC[sim_params.cells] + hWithBC[sim_params.cells + 1];
		etaWest[0] = etaTemp[1] - hWithBC[1] + hWithBC[0];

		for (int i = 0; i < sim_params.cells + 1; i++)
		{
			real uWest = (hWest[i] <= solver_params.tol_dry) ? 0 : qWest[i] / hWest[i];
			real uEast = (hEast[i] <= solver_params.tol_dry) ? 0 : qEast[i] / hEast[i];

			real zWest = etaWest[i] - hWest[i];
			real zEast = etaEast[i] - hEast[i];

			real zIntermediate = max(zWest, zEast);

			deltaWest[i] = max(C(0.0), -(etaWest[i] - zIntermediate));
			deltaEast[i] = max(C(0.0), -(etaEast[i] - zIntermediate));

			hWestStar[i] = max(C(0.0), etaWest[i] - zIntermediate);
			qWestStar[i] = uWest * hWestStar[i];

			hEastStar[i] = max(C(0.0), etaEast[i] - zIntermediate);
			qEastStar[i] = uEast * hEastStar[i];

			zWestStar[i] = zIntermediate - deltaWest[i];
			zEastStar[i] = zIntermediate - deltaEast[i];
		}

		// initialising numerical fluxes
		fluxHLL(sim_params, solver_params, hWestStar, hEastStar, qWestStar, qEastStar, massFlux, momentumFlux);

		for (int i = 0; i < sim_params.cells; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			hBar[i] = (hWestStar[i + 1] + hEastStar[i]) / 2;

			// 1st order projection
			zBar[i] = (zWestStar[i + 1] - zEastStar[i]) / (2 * sqrt(C(3.0)));
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

				hWithBC[i] += dt * massIncrement;
				qWithBC[i] = (hWithBC[i] <= solver_params.tol_dry) ? 0 : qWithBC[i] + dt * momentumIncrement;
			}
		}

		// CFL time step adjustmenttotal_mass 
		dt = 1e9;
		real total_mass = 0;

		for (int i = 1; i < sim_params.cells + 1; i++)
		{
			if (hWithBC[i] > solver_params.tol_dry)
			{
				real u = qWithBC[i] / hWithBC[i];
				real dtCFL = solver_params.CFL * dx / (abs(u) + sqrt(solver_params.g * hWithBC[i]));
				dt = min(dt, dtCFL);
			}

			total_mass  += hWithBC[i] * dx;
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
		test << (xInt[i] + xInt[i + 1]) / 2 << "," << qWithBC[i + 1] << "," << zWithBC[i + 1] << "," << max(zWithBC[i + 1], hWithBC[i + 1] + zWithBC[i + 1]) << endl;
	}

	test.close();

	// delete buffers
	delete[] xInt;

	delete[] qInt;
	delete[] hInt;
	delete[] zInt;

	delete[] qWithBC;
	delete[] hWithBC;
	delete[] zWithBC;

	delete[] dry;

	delete[] etaTemp;

	delete[] qEast;
	delete[] hEast;
	delete[] etaEast;

	delete[] qWest;
	delete[] hWest;
	delete[] etaWest;

	delete[] qWestStar;
	delete[] hWestStar;

	delete[] qEastStar;
	delete[] hEastStar;

	delete[] zWestStar;
	delete[] zEastStar;

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