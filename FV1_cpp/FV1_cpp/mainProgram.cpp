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

	// for-loop index
	int i;

	int steps = 0;

	std::cout << "Please enter a number between 1 and 6 to select a test case.\n"
		         "1: Wet dam break\n"
		         "2: Dry dam break\n"
		         "3: Dry dam break with friction\n"
		         "4: Wet lake-at-rest (C-property)\n"
		         "5: Wet/dry lake-at-rest\n"
		         "6: Building overtopping\n";

	int test_case_selection;

	std::cin >> test_case_selection;

	if (!std::cin || test_case_selection > 6 || test_case_selection < 1)
	{
		std::cout << "Error: please rerun and enter a number between 1 and 6. Exiting program.\n";

		return -1;
	}

	std::cout << "Please enter the number of cells.\n";

	int number_of_cells;

	std::cin >> number_of_cells;

	if (!std::cin || number_of_cells < 1)
	{
		std::cout << "Error: please rerun and enter a integer value. Exiting program.\n";

		return -1;
	}

	SimulationParameters simulationParameters = set_simulation_parameters_for_test_case(test_case_selection, number_of_cells);
	SolverParameters solverParameters = set_solver_parameters_for_test_case();
	BoundaryConditions bcs = set_boundary_conditions_for_test_case(test_case_selection);

	// coarsest cell size
	real dx = (simulationParameters.xmax - simulationParameters.xmin) / simulationParameters.cells;

	// allocate buffer for interfaces
	real* xInt = new real[simulationParameters.cells + 1];

	// initialise baseline mesh
	for (i = 0; i < simulationParameters.cells + 1; i++)
	{
		xInt[i] = simulationParameters.xmin + i * dx;
	}

	// allocate buffers for flow nodes
	real* qInt = new real[simulationParameters.cells + 1];
	real* hInt = new real[simulationParameters.cells + 1];
	real* zInt = new real[simulationParameters.cells + 1];

	// initial interface values
	switch (test_case_selection)
	{
	case 1: case 2: case 3:
		for (i = 0; i < simulationParameters.cells + 1; i++)
		{
			zInt[i] = 0;
			hInt[i] = hInitialOvertopping(bcs, zInt[i], xInt[i]);
			qInt[i] = qInitial(bcs, xInt[i]);
		}
		break;
	case 4: case 5:
		for (i = 0; i < simulationParameters.cells + 1; i++)
		{
			zInt[i] = bedDataConservative(xInt[i]);
			hInt[i] = hInitialConservative(bcs, zInt[i], xInt[i]);
			qInt[i] = qInitial(bcs, xInt[i]);
		}
		break;
	case 6:
		for (i = 0; i < simulationParameters.cells + 1; i++)
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
	real* qWithBC = new real[simulationParameters.cells + 2];
	real* hWithBC = new real[simulationParameters.cells + 2];
	real* zWithBC = new real[simulationParameters.cells + 2];

	// project to find the modes
	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		qWithBC[i] = (qInt[i - 1] + qInt[i]) / 2;
		hWithBC[i] = (hInt[i - 1] + hInt[i]) / 2;
		zWithBC[i] = (zInt[i - 1] + zInt[i]) / 2;
	}

	// allocate true/false buffer for dry cells
	int* dry = new int[simulationParameters.cells + 2];

	real* etaTemp = new real[simulationParameters.cells + 2];

	// allocating buffers for eastern and western interface values
	real* qEast = new real[simulationParameters.cells + 1];
	real* hEast = new real[simulationParameters.cells + 1];
	real* etaEast = new real[simulationParameters.cells + 1];

	real* qWest = new real[simulationParameters.cells + 1];
	real* hWest = new real[simulationParameters.cells + 1];
	real* etaWest = new real[simulationParameters.cells + 1];

	// allocating buffers for positivity preserving nodes
	real* qEastStar = new real[simulationParameters.cells + 1];
	real* hEastStar = new real[simulationParameters.cells + 1];

	real* qWestStar = new real[simulationParameters.cells + 1];
	real* hWestStar = new real[simulationParameters.cells + 1];

	real* zWestStar = new real[simulationParameters.cells + 1];
	real* zEastStar = new real[simulationParameters.cells + 1];
	
	real* deltaWest = new real[simulationParameters.cells + 1];
	real* deltaEast = new real[simulationParameters.cells + 1];

	// allocating buffers for numerical fluxes from HLL solver
	real* massFlux = new real[simulationParameters.cells + 1];
	real* momentumFlux = new real[simulationParameters.cells + 1];

	// allocating buffers for positivity preserving MODES
	real* hBar = new real[simulationParameters.cells];
	real* zBar = new real[simulationParameters.cells];

	real timeNow = 0;
	real dt = C(1e-3);

	int firstTimeStep = 1;

	// file stream for recording cumulative clock time vs sim time
	ofstream data;

	data.open("clock_time_vs_sim_time.csv");

	data.precision(15);

	data << "sim_time,clock_time" << std::endl;

	while (timeNow < simulationParameters.simulationTime)
	{
		timeNow += dt;
		
		if (timeNow - simulationParameters.simulationTime > 0)
		{
			timeNow -= dt;
			dt = simulationParameters.simulationTime - timeNow;
			timeNow += dt;
		}

		// adding ghost boundary conditions
		qWithBC[0] = (bcs.qxImposedUp > 0) ? bcs.qxImposedUp : bcs.reflectUp * qWithBC[1];
		qWithBC[simulationParameters.cells + 1] = (bcs.qxImposedDown > 0) ? bcs.qxImposedDown : bcs.reflectDown * qWithBC[simulationParameters.cells];

		hWithBC[0] = (bcs.hImposedUp > 0) ? bcs.hImposedUp : hWithBC[1];
		hWithBC[simulationParameters.cells + 1] = (bcs.hImposedDown > 0) ? bcs.hImposedDown : hWithBC[simulationParameters.cells];

		zWithBC[0] = zWithBC[1];
		zWithBC[simulationParameters.cells + 1] = zWithBC[simulationParameters.cells];

		if (simulationParameters.manning > 0)
		{
			for (i = 1; i < simulationParameters.cells + 1; i++)
			{
				qWithBC[i] += frictionImplicit(simulationParameters, solverParameters, dt, hWithBC[i], qWithBC[i]);
			}
		}

		// ghost cells are always 0 i.e. false/wet
		dry[0] = false;
		dry[simulationParameters.cells + 1] = false;

		// initialising dry vs wet cells, ignore ghost cells
		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			real hLocal = hWithBC[i];
			real hBackward = hWithBC[i - 1];
			real hForward = hWithBC[i + 1];

			real hMax = max(hLocal, hBackward);
			hMax = max(hForward, hMax);

			dry[i] = (hMax <= solverParameters.tolDry);
		}

		for (i = 0; i < simulationParameters.cells + 2; i++)
		{
			etaTemp[i] = hWithBC[i] + zWithBC[i];
		}

		// initialising interface values
		for (i = 0; i < simulationParameters.cells + 1; i++)
		{
			qEast[i] = qWithBC[i + 1];
			hEast[i] = hWithBC[i + 1];
			etaEast[i] = etaTemp[i + 1];

			qWest[i] = qWithBC[i];
			hWest[i] = hWithBC[i];
			etaWest[i] = etaTemp[i];
		}

		// correcting downwind and upwind eta values
		etaEast[simulationParameters.cells] = etaTemp[simulationParameters.cells] - hWithBC[simulationParameters.cells] + hWithBC[simulationParameters.cells + 1];
		etaWest[0] = etaTemp[1] - hWithBC[1] + hWithBC[0];

		for (int i = 0; i < simulationParameters.cells + 1; i++)
		{
			real uWest = (hWest[i] <= solverParameters.tolDry) ? 0 : qWest[i] / hWest[i];
			real uEast = (hEast[i] <= solverParameters.tolDry) ? 0 : qEast[i] / hEast[i];

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
		fluxHLL(simulationParameters, solverParameters, hWestStar, hEastStar, qWestStar, qEastStar, massFlux, momentumFlux);

		for (int i = 0; i < simulationParameters.cells; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			hBar[i] = (hWestStar[i + 1] + hEastStar[i]) / 2;

			// 1st order projection
			zBar[i] = (zWestStar[i + 1] - zEastStar[i]) / (2 * sqrt(C(3.0)));
		}

		// FV1 operator increment, skip ghosts cells
		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			// skip increment in dry cells
			if (!dry[i])
			{
				real massIncrement = -(1 / dx) * (massFlux[i] - massFlux[i - 1]);
				real momentumIncrement = -(1 / dx) * (momentumFlux[i] - momentumFlux[i - 1] + 2 * sqrt(C(3.0)) * solverParameters.g * hBar[i - 1] * zBar[i - 1]);

				real a = momentumFlux[i] - momentumFlux[i - 1];
				real b = 2 * sqrt(C(3.0)) * solverParameters.g * hBar[i - 1] * zBar[i - 1];

				hWithBC[i] += dt * massIncrement;
				qWithBC[i] = (hWithBC[i] <= solverParameters.tolDry) ? 0 : qWithBC[i] + dt * momentumIncrement;
			}
		}

		// CFL time step adjustmenttotal_mass 
		dt = 1e9;
		real total_mass = 0;

		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			if (hWithBC[i] > solverParameters.tolDry)
			{
				real u = qWithBC[i] / hWithBC[i];
				real dtCFL = solverParameters.CFL * dx / (abs(u) + sqrt(solverParameters.g * hWithBC[i]));
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

	for (i = 0; i < simulationParameters.cells; i++)
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