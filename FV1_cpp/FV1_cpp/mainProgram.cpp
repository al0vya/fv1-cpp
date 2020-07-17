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
#include "bedDataConservative.h"
#include "qInitial.h"
#include "hInitial.h"
#include "frictionImplicit.h"
#include "uFaceValues.h"
#include "fluxHLL.h"
#include "fluxHLLNew.h"

using namespace std;

int main()
{
	clock_t start = clock();

	// quintessential for-loop index
	int i;

	int step = 0;

	SimulationParameters simulationParameters;
	simulationParameters.cells = 512;
	simulationParameters.xmin = 0;
	simulationParameters.xmax = 50;
	simulationParameters.simulationTime = C(100.0);
	simulationParameters.manning = C(0.0);

	SolverParameters solverParameters;
	solverParameters.CFL = C(0.33);
	solverParameters.tolDry = C(1e-3);
	solverParameters.g = C(9.80665);

	BoundaryConditions bcs;
	bcs.hl = C(2.0);
	bcs.hr = C(2.0);

	bcs.ql = C(0.0);
	bcs.qr = C(0.0);

	bcs.reflectUp = C(1.0);
	bcs.reflectDown = C(1.0);

	bcs.hImposedUp = C(0.0);
	bcs.qxImposedUp = C(0.0);

	bcs.hImposedDown = C(0.0);
	bcs.qxImposedDown = C(0.0);

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
	for (i = 0; i < simulationParameters.cells + 1; i++)
	{
		zInt[i] = bedDataConservative(xInt[i]);
		hInt[i] = hInitial(bcs, zInt[i], xInt[i]);
		qInt[i] = qInitial(bcs, xInt[i]);
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

		// extract upwind and downwind modes
		real hWestUpwind = hWithBC[0];
		real qWestUpwind = qWithBC[0];

		real hEastDownwind = hWithBC[simulationParameters.cells + 1];
		real qEastDownwind = qWithBC[simulationParameters.cells + 1];

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
		etaEast[simulationParameters.cells] = etaTemp[simulationParameters.cells] - hWithBC[simulationParameters.cells] + hEastDownwind;
		etaWest[0] = etaTemp[1] - hWithBC[1] + hWestUpwind;

		for (int i = 0; i < simulationParameters.cells + 1; i++)
		{
			real zWest = etaWest[i] - hWest[i];
			real zEast = etaEast[i] - hEast[i];
			
			real uWest = (hWest[i] <= solverParameters.tolDry) ? 0 : qWest[i] / hWest[i];
			real uEast = (hEast[i] <= solverParameters.tolDry) ? 0 : qEast[i] / hEast[i];

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

		if (firstTimeStep)
		{
			ofstream arrays;

			arrays.open("arrays.txt");

			arrays << "xInt, qInt, hInt, zInt, qWest, hWest, etaWest, qEast, hEast, etaEast, "
				"qWestStar, hWestStar, zWestStar, qEastStar, hEastStar, zEastStar, deltaWest, deltaEast, massFlux, momentumFlux" << endl;

			for (int i = 0; i < simulationParameters.cells + 1; i++)
			{
				arrays << xInt[i] << "," << qInt[i] << "," << hInt[i] << "," << zInt[i] << "," 
					<< qWest[i] << "," << hWest[i] << "," << etaWest[i] << "," << qEast[i] << "," << hEast[i] << "," << etaEast[i] << ","
					<< qWestStar[i] << "," << hWestStar[i] << "," << zWestStar[i] << "," << qEastStar[i] << "," << hEastStar[i] << "," 
					<< zEastStar[i] << "," << deltaWest[i] << "," << deltaEast[i] << ","
					<< massFlux[i] << "," << momentumFlux[i] << endl;
			}

			arrays << endl;

			arrays << "qWithBC, hWithBC, zWithBC, etaTemp, dry" << endl;

			for (i = 0; i < simulationParameters.cells + 2; i++)
			{
				arrays << qWithBC[i] << "," << hWithBC[i] << "," << zWithBC[i] << "," << etaTemp[i] << "," << dry[i] << endl;
			}

			arrays << endl;

			arrays << "hBar, zBar" << endl;

			for (i = 0; i < simulationParameters.cells; i++)
			{
				arrays << hBar[i] << "," << zBar[i] << endl;
			}

			firstTimeStep = 0;
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
				qWithBC[i] = hWithBC[i] <= solverParameters.tolDry ? 0 : qWithBC[i] + dt * momentumIncrement;
			}
		}

		// CFL time step adjustment
		dt = 1e9;
		real totalMass = 0;

		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			if (hWithBC[i] > solverParameters.tolDry)
			{
				real u = qWithBC[i] / hWithBC[i];
				real dtCFL = solverParameters.CFL * dx / (abs(u) + sqrt(solverParameters.g * hWithBC[i]));
				dt = min(dt, dtCFL);
			}

			totalMass += hWithBC[i] * dx;
		}

		//step++;

		printf("Mass: %.17g, time: %f s\n", totalMass, timeNow);

		
	}

	ofstream test;

	test.open("FV1_data.txt");

	test << "x,q,x,z,eta" << endl;

	for (i = 0; i < simulationParameters.cells; i++)
	{
		test << (xInt[i] + xInt[i + 1]) / 2 << "," << qWithBC[i + 1] << "," << (xInt[i] + xInt[i + 1]) / 2 << "," << zWithBC[i + 1] << "," << hWithBC[i + 1] + zWithBC[i + 1] << endl;
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

	clock_t end = clock();

	real time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", time);

	return 0;
}