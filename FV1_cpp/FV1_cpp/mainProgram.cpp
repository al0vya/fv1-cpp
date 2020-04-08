#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>

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

using namespace std;

int main()
{
	// quintessential for-loop index
	int i;

	int step = 0;

	SimulationParameters simulationParameters;
	simulationParameters.cells = 100;
	simulationParameters.xmin = 0;
	simulationParameters.xmax = 50;
	simulationParameters.simulationTime = C(20.);
	simulationParameters.manning = C(0.02);

	SolverParameters solverParameters;
	solverParameters.CFL = C(0.33);
	solverParameters.tolDry = C(1e-3);
	solverParameters.g = C(9.80665);

	BoundaryConditions bcs;
	bcs.hl = C(6.0);
	bcs.hr = C(0.0);

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
	real* x_int = new real[simulationParameters.cells + 1];

	// initialise baseline mesh
	for (i = 0; i < simulationParameters.cells + 1; i++)
	{
		x_int[i] = simulationParameters.xmin + i * dx;
	}

	// allocate buffers for flow nodes
	real* q_int = new real[simulationParameters.cells + 1];
	real* h_int = new real[simulationParameters.cells + 1];
	real* z_int = new real[simulationParameters.cells + 1];

	// initial interface values
	for (i = 0; i < simulationParameters.cells + 1; i++)
	{
		z_int[i] = bedDataConservative(x_int[i]);
		h_int[i] = hInitial(bcs, z_int[i], x_int[i]);
		q_int[i] = qInitial(bcs, x_int[i]);
	}

	// allocate buffers for flow modes with ghost BCs
	real* qWithBC = new real[simulationParameters.cells + 2];
	real* hWithBC = new real[simulationParameters.cells + 2];
	real* zWithBC = new real[simulationParameters.cells + 2];

	// project to find the modes
	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		qWithBC[i] = (q_int[i - 1] + q_int[i]) / 2;
		hWithBC[i] = (h_int[i - 1] + h_int[i]) / 2;
		zWithBC[i] = (z_int[i - 1] + z_int[i]) / 2;
	}

	// allocate true/false buffer for dry cells
	bool* dry = new bool[simulationParameters.cells + 2];

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

	real* zStarIntermediate = new real[simulationParameters.cells + 1];
	real* zStar = new real[simulationParameters.cells + 1];
	 
	real* uWest = new real[simulationParameters.cells + 1];
	real* uEast = new real[simulationParameters.cells + 1];

	real* delta = new real[simulationParameters.cells + 1];

	// allocating buffers for numerical fluxes from HLL solver
	real* massFlux = new real[simulationParameters.cells + 1];
	real* momentumFlux = new real[simulationParameters.cells + 1];

	// allocating buffers for positivity preserving MODES
	real* hBar = new real[simulationParameters.cells];
	real* zBar = new real[simulationParameters.cells];

	real timeNow = 0;
	real dt = C(1e-4);

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
		real qUp = bcs.reflectUp * qWithBC[1];
		if (bcs.qxImposedUp > 0)
		{
			qUp = bcs.qxImposedUp;
		}

		real qDown = bcs.reflectDown * qWithBC[simulationParameters.cells];
		if (bcs.qxImposedDown > 0)
		{
			qDown = bcs.qxImposedDown;
		}

		qWithBC[0] = qUp;
		qWithBC[simulationParameters.cells + 1] = qDown; // recall there are cells + 2 elements inc BCs

		real hUp = hWithBC[1];
		if (bcs.hImposedUp > 0)
		{
			hUp = bcs.hImposedUp;
		}

		real hDown = hWithBC[simulationParameters.cells];
		if (bcs.hImposedDown > 0)
		{
			hDown = bcs.hImposedDown;
		}

		hWithBC[0] = hUp;
		hWithBC[simulationParameters.cells + 1] = hDown;

		real zUp = zWithBC[1];
		real zDown = zWithBC[simulationParameters.cells];

		zWithBC[0] = zUp;
		zWithBC[simulationParameters.cells + 1] = zDown;

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

			// dry[] hasn't been initialised so else statement also needed
			if (hMax <= solverParameters.tolDry)
			{
				dry[i] = true;
			}
			else
			{
				dry[i] = false;
			}
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
			// initialising velocity interface values
			uWest[i] = uFaceValues(solverParameters, qWest[i], hWest[i]);
			uEast[i] = uFaceValues(solverParameters, qEast[i], hEast[i]);

			// intermediate calculations
			real a = etaWest[i] - hWest[i];
			real b = etaEast[i] - hEast[i];
			zStarIntermediate[i] = max(a, b);
			a = etaWest[i] - zStarIntermediate[i];
			b = etaEast[i] - zStarIntermediate[i];

			// positivity-preserving nodes
			hWestStar[i] = max(C(0.0), a);
			hEastStar[i] = max(C(0.0), b);

			delta[i] = max(C(0.0), -a) + max(C(0.0), -b);

			qWestStar[i] = uWest[i] * hWestStar[i];
			qEastStar[i] = uEast[i] * hEastStar[i];

			zStar[i] = zStarIntermediate[i] - delta[i];
		}

		// initialising numerical fluxes
		fluxHLL(simulationParameters, solverParameters, hWestStar, hEastStar, qWestStar, qEastStar, uWest, uEast, massFlux, momentumFlux);

		for (int i = 0; i < simulationParameters.cells; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			hBar[i] = (hEastStar[i] + hWestStar[i + 1]) / 2;

			// first order projection
			zBar[i] = (zStar[i + 1] - zStar[i]) / (2 * sqrt(C(3.0)));
		}

		// FV1 operator increment, skip ghosts cells
		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			// skip increment in dry cells
			if (dry[i])
			{
				continue;
			}
			else
			{
				real massIncrement = -(1 / dx) * (massFlux[i] - massFlux[i - 1]);
				real momentumIncrement = -(1 / dx) * (momentumFlux[i] - momentumFlux[i - 1] + 2 * sqrt(C(3.0)) * solverParameters.g * hBar[i - 1] * zBar[i - 1]);

				hWithBC[i] += dt * massIncrement;
				qWithBC[i] += dt * momentumIncrement;
			}

			if (hWithBC[i] <= solverParameters.tolDry)
			{
				qWithBC[i] = 0;
			}
		}

		// CFL time step adjustment
		dt = 1e9;

		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			if (hWithBC[i] <= solverParameters.tolDry)
			{
				continue;
			}
			else
			{
				real u = qWithBC[i] / hWithBC[i];
				real dtCFL = solverParameters.CFL * dx / (abs(u) + sqrt(solverParameters.g * hWithBC[i]));
				dt = min(dt, dtCFL);
			}
		}

		step++;
		
		//printf("%f s\n", timeNow);
	}

	for (i = 1; i < simulationParameters.cells + 1; i++)
	{
		printf("%0.2f\n", hWithBC[i] + zWithBC[i]);
	}

	printf("\n");

	// delete buffers
	delete[] x_int;

	delete[] q_int;
	delete[] h_int;
	delete[] z_int;

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

	delete[] zStarIntermediate;
	delete[] zStar;

	delete[] uWest;
	delete[] uEast;

	delete[] delta;

	delete[] massFlux;
	delete[] momentumFlux;

	delete[] hBar;
	delete[] zBar;

	return 0;
}