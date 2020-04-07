#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include "baselineMesh.h"
#include "real.h"
#include "simulationParameters.h"

using namespace std;

typedef struct BoundaryConditions
{
	real hl;
	real hr;

	real ql;
	real qr;
	
	real reflectUp;
	real reflectDown;

	real hImposedUp;
	real qxImposedUp;

	real hImposedDown;
	real qxImposedDown;
} BoundaryConditions;

/*
requires /std:c++latest
the (slight) advantage is bcs will be initialised upon construction
there is never a time when bcs is (partially) uninitialised

BoundaryConditions bcs = {
	.hl = C(2.0),
	.hr = C(0.0),
	.ql = C(0.0),
	.qr = C(0.0),
	.reflectUp = C(0.0),
	.reflectDown = C(0.0),
	.hImposedUp = C(0.0),
	.qxImposedUp = C(0.0),
	.hImposedDown = C(0.0),
	.qxImposedDown = C(0.0)
};
*/

typedef struct SolverParameters
{
	real CFL;
	real tolDry;
	real g;
} SolverParameters;

// declaring helper functions
void bedDataConservative(SimulationParameters simulationParameters, real* x_int, real* z_int);
void bedDataDamBreak(SimulationParameters simulationParameters, real* x_int, real* z_int);
void qInitialDamBreak(SimulationParameters simulationParameters, BoundaryConditions bcs, real* x_int, real* q_int);
void hInitialDamBreak(SimulationParameters simulationParameters, BoundaryConditions bcs, real* z_int, real* x_int, real* h_int);
void modalProjectionZeroOrder(SimulationParameters simulationParameters, real* u_int, real* u);
void modalProjectionFirstOrder(SimulationParameters simulationParameters, real* u_int, real* u);
void qAddGhostBoundaryConditions(SimulationParameters simulationParameters, BoundaryConditions bcs, real* q, real* qWithBC);
void hAddGhostBoundaryConditions(SimulationParameters simulationParameters, BoundaryConditions bcs, real* h, real* hWithBC);
void zAddGhostBoundaryConditions(SimulationParameters simulationParameters, real* z, real* zWithBC);
void hWestFaceValues(SimulationParameters simulationParameters, real hWestUpwind, real* hTemp, real* hWest);
void qWestFaceValues(SimulationParameters simulationParameters, real qWestUpwind, real* qTemp, real* qWest);
void etaWestFaceValues(SimulationParameters simulationParameters, real hWestUpwind, real* hTemp, real* etaTemp, real* etaWest);
void hEastFaceValues(SimulationParameters simulationParameters, real hEastDownwind, real* hTemp, real* hEast);
void qEastFaceValues(SimulationParameters simulationParameters, real qEastDownwind, real* qTemp, real* qEast);
void etaEastFaceValues(SimulationParameters simulationParameters, real hEastDownwind, real* hTemp, real* etaTemp, real* etaEast);
void uFaceValues(SimulationParameters simulationParameters, SolverParameters solverParameters, real* qFace, real* hFace, real* uFace);
void zStarIntermediateValues(SimulationParameters simulationParameters, real* etaWest, real* etaEast, real* hWest, real* hEast, real* zStarIntermediate);
void deltaValues(SimulationParameters simulationParameters, real* etaFace, real* zStarIntermediate, real* deltaFace);
void hStarValues(SimulationParameters simulationParameters, real* etaFaceValue, real* zStarIntermediate, real* hStar);
void qStarValues(SimulationParameters simulationParameters, real* uFaceValue, real* hStar, real* qStar);
void zStarValues(SimulationParameters simulationParameters, real* deltaWest, real* deltaEast, real* zStarIntermediate, real* zStar);
void fluxHLL(SimulationParameters simulationParameters, SolverParameters solverParameters, real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* uWest, real* uEast, real* massFlux, real* momentumFlux);
void hBarValues(SimulationParameters simulationParameters, real* hWestStar, real* hEastStar, real* hBar);
void momentumFV1OperatorValues(SimulationParameters simulationParameters, SolverParameters solverParameters, real dx, real* hBar, real* zBar, real* momentumFlux, real* momentumFV1Operator);
void frictionImplicit(SimulationParameters simulationParameters, SolverParameters solverParameters, real dt, real* hWithBC, real* qWithBC);

int main()
{
	// quintessential for-loop index
	int i;

	int step = 0;

	SimulationParameters simulationParameters;
	simulationParameters.cells = 10;
	simulationParameters.xmin = 0;
	simulationParameters.xmax = 50;
	simulationParameters.simulationTime = C(500.0);
	simulationParameters.manning = C(0.02);

	SolverParameters solverParameters;
	solverParameters.CFL = C(0.33);
	solverParameters.tolDry = C(1e-3);
	solverParameters.g = C(9.80665);

	real u, dtCFL;

	// coarsest cell size
	real dx = (simulationParameters.xmax - simulationParameters.xmin) / simulationParameters.cells;

	// allocate buffers for the cell centres and interfaces
	real* x = new real[simulationParameters.cells];
	real* x_int = new real[simulationParameters.cells + 1];

	// initialise baseline mesh
	baselineMesh(simulationParameters, dx, x, x_int);

	// allocate buffers for flow modes q, h, z
	real* q = new real[simulationParameters.cells];
	real* h = new real[simulationParameters.cells];
	real* z = new real[simulationParameters.cells];

	// allocate buffers for the corresponding interface Valuess
	real* q_int = new real[simulationParameters.cells + 1];
	real* h_int = new real[simulationParameters.cells + 1];
	real* z_int = new real[simulationParameters.cells + 1];

	BoundaryConditions bcs;
	bcs.hl = C(12.0);
	bcs.hr = C(0.0);

	bcs.ql = C(0.0);
	bcs.qr = C(0.0);

	bcs.reflectUp = C(1.0);
	bcs.reflectDown = C(-1.0);

	bcs.hImposedUp = C(0.0);
	bcs.qxImposedUp = C(0.0);

	bcs.hImposedDown = C(0.0);
	bcs.qxImposedDown = C(0.0);

	// initialise the interface Valuess
	bedDataConservative(simulationParameters, x_int, z_int);
	qInitialDamBreak(simulationParameters, bcs, x_int, q_int);
	hInitialDamBreak(simulationParameters, bcs, z_int, x_int, h_int);

	// project to find the modes
	modalProjectionZeroOrder(simulationParameters, q_int, q);
	modalProjectionZeroOrder(simulationParameters, h_int, h);
	modalProjectionZeroOrder(simulationParameters, z_int, z);

	// allocate buffers for modes with BCs
	real* qWithBC = new real[simulationParameters.cells + 2];
	real* hWithBC = new real[simulationParameters.cells + 2];
	real* zWithBC = new real[simulationParameters.cells + 2];

	// allocate buffers for RK2 step
	real* qWithBCNew = new real[simulationParameters.cells + 2];
	real* hWithBCNew = new real[simulationParameters.cells + 2];

	real* qTemp = new real[simulationParameters.cells + 2];
	real* hTemp = new real[simulationParameters.cells + 2];
	real* zTemp = new real[simulationParameters.cells + 2];

	// placeholder variables for substitution
	real a, b, c;

	// allocate true/false buffer for dry cells
	bool* dry = new bool[simulationParameters.cells + 2];

	// placeholder variables for checking
	real hLocal, hBackward, hForward, hMax;

	real* etaTemp = new real[simulationParameters.cells + 2];

	// allocating buffers for eastern and western interface Valuess
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

	real* deltaWest = new real[simulationParameters.cells + 1];
	real* deltaEast = new real[simulationParameters.cells + 1];

	// allocating buffers for numerical fluxes with HLL solver
	real* massFlux = new real[simulationParameters.cells + 1];
	real* momentumFlux = new real[simulationParameters.cells + 1];

	// allocating buffers for positivity preserving MODES
	real* hBar = new real[simulationParameters.cells];
	real* zBar = new real[simulationParameters.cells];

	// allocating buffer for fv1Operator values
	real* massFV1Operator = new real[simulationParameters.cells];
	real* momentumFV1Operator = new real[simulationParameters.cells];

	// WHILE LOOP STARTS FROM HERE //
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

		// initialise modes with BCs
		qAddGhostBoundaryConditions(simulationParameters, bcs, q, qWithBC);
		hAddGhostBoundaryConditions(simulationParameters, bcs, h, hWithBC);
		zAddGhostBoundaryConditions(simulationParameters, z, zWithBC);

		// extract upwind and downwind modes
		real hWestUpwind = hWithBC[0];
		real qWestUpwind = qWithBC[0];

		real hEastDownwind = hWithBC[simulationParameters.cells + 1];
		real qEastDownwind = qWithBC[simulationParameters.cells + 1];

		// initialise RK2 buffers
		for (i = 0; i < simulationParameters.cells + 2; i++)
		{
			a = qWithBC[i];
			b = hWithBC[i];
			c = zWithBC[i];

			qWithBCNew[i] = a;
			hWithBCNew[i] = b;

			qTemp[i] = a;
			hTemp[i] = b;
			zTemp[i] = c;
		}

		if (simulationParameters.manning > 0)
		{
			frictionImplicit(simulationParameters, solverParameters, dt, hWithBC, qWithBC);
		}

		// ghost cells are always 0 i.e. false/wet
		dry[0] = false;
		dry[simulationParameters.cells + 1] = false;

		// initialising dry vs wet cells, ignore ghost cells
		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			hLocal = hTemp[i];
			hBackward = hTemp[i - 1];
			hForward = hTemp[i + 1];

			hMax = max(hLocal, hBackward);
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
			etaTemp[i] = hTemp[i] + zTemp[i];
		}
		/*
		__global__ fvoperator()
		{
		
		for ( int i = 0; i < cells; i++) // in CUDA, this would be a grid-stride loop running in parallel (each thread picks a different i value)
		{
			// do loads of stuff that's local to this cell
			// this isn't entirely possible: there's fluxes of mass and momentum between cells
			// so there has to be some synchronisation
		}
		*/

		// initialising interface Valuess
		qEastFaceValues(simulationParameters, qEastDownwind, qTemp, qEast);
		hEastFaceValues(simulationParameters, hEastDownwind, hTemp, hEast);
		etaEastFaceValues(simulationParameters, hEastDownwind, hTemp, etaTemp, etaEast);

		qWestFaceValues(simulationParameters, qWestUpwind, qTemp, qWest);
		hWestFaceValues(simulationParameters, hWestUpwind, hTemp, hWest);
		etaWestFaceValues(simulationParameters, hWestUpwind, hTemp, etaTemp, etaWest);

		// initialising velocity interface values
		uFaceValues(simulationParameters, solverParameters, qWest, hWest, uWest);
		uFaceValues(simulationParameters, solverParameters, qEast, hEast, uEast);

		zStarIntermediateValues(simulationParameters, etaWest, etaEast, hWest, hEast, zStarIntermediate);

		deltaValues(simulationParameters, etaWest, zStarIntermediate, deltaWest);
		deltaValues(simulationParameters, etaEast, zStarIntermediate, deltaEast);

		// initialising positivity preserving nodes
		hStarValues(simulationParameters, etaWest, zStarIntermediate, hWestStar);
		hStarValues(simulationParameters, etaEast, zStarIntermediate, hEastStar);

		qStarValues(simulationParameters, uWest, hWestStar, qWestStar);
		qStarValues(simulationParameters, uEast, hEastStar, qEastStar);

		zStarValues(simulationParameters, deltaWest, deltaEast, zStarIntermediate, zStar);

		// initialising numerical fluxes
		fluxHLL(simulationParameters, solverParameters, hWestStar, hEastStar, qWestStar, qEastStar, uWest, uEast, massFlux, momentumFlux);

		hBarValues(simulationParameters, hWestStar, hEastStar, hBar);

		modalProjectionFirstOrder(simulationParameters, zStar, zBar);

		momentumFV1OperatorValues(simulationParameters, solverParameters, dx, hBar, zBar, momentumFlux, momentumFV1Operator);

		// FV1 operator increment, skip ghosts cells
		for (i = 1; i < simulationParameters.cells + 1; i++)
		{
			// skip increment in dry cells
			if (dry[i])
			{
				hWithBCNew[i] = hWithBC[i];
				qWithBCNew[i] = qWithBC[i];
			}
			else
			{
				real massIncrement = -(1 / dx) * (massFlux[i] - massFlux[i - 1]);		

				hWithBCNew[i] = hWithBC[i] + dt * massIncrement;
				qWithBCNew[i] = qWithBC[i] + dt * momentumFV1Operator[i - 1]; // TODO: do the same thing
			}

			if (hWithBCNew[i] <= solverParameters.tolDry)
			{
				qWithBCNew[i] = 0;
			}
		}

		// update and CFL time step adjustment
		dt = 1e9;

		for (i = 0; i < simulationParameters.cells; i++)
		{
			h[i] = hWithBCNew[i + 1];
			q[i] = qWithBCNew[i + 1];

			if (h[i] <= solverParameters.tolDry)
			{
				continue;
			}
			else
			{
				u = q[i] / h[i];
				dtCFL = solverParameters.CFL * dx / (abs(u) + std::sqrt(solverParameters.g * h[i]));
				dt = min(dt, dtCFL);
			}
		}

		step++;

		for (i = 0; i < simulationParameters.cells; i++)
		{
			printf("%f, ", h[i] + z[i]);
		}
		printf("%f s\n", timeNow);
	}

	printf("Finished.\n");

	// delete buffers
	delete[] x;
	delete[] x_int;

	delete[] q;
	delete[] h;
	delete[] z;

	delete[] q_int;
	delete[] h_int;
	delete[] z_int;

	delete[] qWithBC;
	delete[] hWithBC;
	delete[] zWithBC;

	delete[] qWithBCNew;
	delete[] hWithBCNew;

	delete[] qTemp;
	delete[] hTemp;
	delete[] zTemp;

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

	delete[] deltaWest;
	delete[] deltaEast;

	delete[] massFlux;
	delete[] momentumFlux;

	delete[] hBar;
	delete[] zBar;

	delete[] massFV1Operator;
	delete[] momentumFV1Operator;

	return 0;
}

// Helper function definitions //

void baselineMesh(SimulationParameters simulationParameters, real dx, real* x, real* x_int)
{
	real a, b;

	x_int[0] = simulationParameters.xmin;
	x_int[simulationParameters.cells] = simulationParameters.xmax;

	for (int i = 1; i < simulationParameters.cells; i++) {
		a = x_int[i - 1];
		b = a + dx;

		x_int[i] = b;
		x[i - 1] = (a + b) / 2;
	}

	// final cell centre since it's not done in the for loop
	x[simulationParameters.cells - 1] = (x_int[simulationParameters.cells] + x_int[simulationParameters.cells - 1]) / 2;
}

void bedDataConservative(SimulationParameters simulationParameters, real* x_int, real* z_int)
{
	real a;

	// for (cells+1) interfaces
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		a = x_int[i];

		if (a >= 22 && a < 25)
		{
			z_int[i] = C(0.05) * a - C(1.1);
		}
		else if (a >= 25 && a <= 28)
		{
			z_int[i] = C(-0.05) * a + C(1.4);
		}
		else if (a > 8 && a < 12)
		{
			z_int[i] = C(0.2) - C(0.05) * pow(a - 10, 2);
		}
		else if (a > 39 && a < 46.5)
		{
			z_int[i] = C(0.3); // whereas here, 0.3 is a double literal, dangerous to cast it to a float
		}
		else
		{
			z_int[i] = 0; // this is safe because you're casting an int literal to a real
		}

		z_int[i] *= 10;
	}
}

// note this function works on the z NODES i.e. interfaces
void bedDataDamBreak(SimulationParameters simulationParameters, real* x_int, real* z_int)
{
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (x_int[i] <= 25)
		{
			z_int[i] = 0;
		}
		else
		{
			z_int[i] = 0;
		}
	}
}

// note this function works on the interfaces
void qInitialDamBreak(SimulationParameters simulationParameters, BoundaryConditions bcs, real* x_int, real* q_int)
{
	// recall for loop is for (cells+1) interfaces
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (x_int[i] <= 32.5)
		{
			q_int[i] = bcs.ql;
		}
		else
		{
			q_int[i] = bcs.qr;
		}
	}
}

void hInitialDamBreak(SimulationParameters simulationParameters, BoundaryConditions bcs, real* z_int, real* x_int, real* h_int)
{
	real etaLeft = bcs.hl;
	real etaRight = bcs.hr;

	// recall for loop is for (cells+1) interfaces
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (x_int[i] <= 25)
		{
			h_int[i] = etaLeft - z_int[i];

			if (h_int[i] < 0)
			{
				h_int[i] = 0;
			}
		}
		else
		{
			h_int[i] = etaRight - z_int[i];

			if (h_int[i] < 0)
			{
				h_int[i] = 0;
			}
		}
	}
}

void modalProjectionZeroOrder(SimulationParameters simulationParameters, real* u_int, real* u)
{
	// for modes, so loops overs cells rather than (cells+1) interfaces
	for (int i = 0; i < simulationParameters.cells; i++)
	{
		u[i] = (u_int[i] + u_int[i + 1]) / 2;
	}
}

void modalProjectionFirstOrder(SimulationParameters simulationParameters, real* u_int, real* u)
{
	for (int i = 0; i < simulationParameters.cells; i++)
	{
		u[i] = (u_int[i + 1] - u_int[i]) / (2 * sqrt(C(3.0)));
	}
}

void  qAddGhostBoundaryConditions(SimulationParameters simulationParameters, BoundaryConditions bcs, real* q, real* qWithBC)
{
	real qUp = bcs.reflectUp * q[0];

	if (bcs.qxImposedUp > 0)
	{
		qUp = bcs.qxImposedUp;
	}

	real qDown = bcs.reflectDown * q[simulationParameters.cells - 1];

	if (bcs.qxImposedDown > 0)
	{
		qDown = bcs.qxImposedDown;
	}

	qWithBC[0] = qUp;

	// there are cells + 2 elements inc BCs
	qWithBC[simulationParameters.cells + 1] = qDown;

	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		qWithBC[i] = q[i - 1];
	}
}

void hAddGhostBoundaryConditions(SimulationParameters simulationParameters, BoundaryConditions bcs, real* h, real* hWithBC)
{
	real hUp = h[0];

	if (bcs.hImposedUp > 0)
	{
		hUp = bcs.hImposedUp;
	}

	real hDown = h[simulationParameters.cells - 1];

	if (bcs.hImposedDown > 0)
	{
		hDown = bcs.hImposedDown;
	}

	hWithBC[0] = hUp;

	// there are cells + 2 elements inc BCs
	hWithBC[simulationParameters.cells + 1] = hDown;

	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		hWithBC[i] = h[i - 1];
	}
}

void zAddGhostBoundaryConditions(SimulationParameters simulationParameters, real* z, real* zWithBC)
{
	real zUp = z[0];
	real zDown = z[simulationParameters.cells - 1];

	zWithBC[0] = zUp;

	// there are cells + 2 elements inc BCs
	zWithBC[simulationParameters.cells + 1] = zDown;

	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		zWithBC[i] = z[i - 1];
	}
}

void hWestFaceValues(SimulationParameters simulationParameters, real hWestUpwind, real* hTemp, real* hWest)
{
	hWest[0] = hWestUpwind;

	// start from 1 on hTemp to avoid ghost cell
	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		hWest[i] = hTemp[i];
	}
}

void qWestFaceValues(SimulationParameters simulationParameters, real qWestUpwind, real* qTemp, real* qWest)
{
	qWest[0] = qWestUpwind;

	// start from 1 on qTemp to avoid ghost cell
	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		qWest[i] = qTemp[i];
	}
}

void etaWestFaceValues(SimulationParameters simulationParameters, real hWestUpwind, real* hTemp, real* etaTemp, real* etaWest)
{
	// indexing shifted by one, since uTemp[0] is the ghost mode
	etaWest[0] = etaTemp[1] - hTemp[1] + hWestUpwind;

	// similarly shifted, start from 1 on temp to avoid ghost cell
	for (int i = 1; i < simulationParameters.cells + 1; i++)
	{
		etaWest[i] = etaTemp[i];
	}
}

void hEastFaceValues(SimulationParameters simulationParameters, real hEastDownwind, real* hTemp, real* hEast)
{
	// [i+1] on hTemp to skip ghost cell
	for (int i = 0; i < simulationParameters.cells; i++)
	{
		hEast[i] = hTemp[i + 1];
	}

	hEast[simulationParameters.cells] = hEastDownwind;
}

void qEastFaceValues(SimulationParameters simulationParameters, real qEastDownwind, real* qTemp, real* qEast)
{
	// [i+1] on qTemp to skip ghost cell
	for (int i = 0; i < simulationParameters.cells; i++)
	{
		qEast[i] = qTemp[i + 1];
	}

	qEast[simulationParameters.cells] = qEastDownwind;
}

void etaEastFaceValues(SimulationParameters simulationParameters, real hEastDownwind, real* hTemp, real* etaTemp, real* etaEast)
{
	for (int i = 0; i < simulationParameters.cells; i++)
	{
		etaEast[i] = etaTemp[i + 1];
	}

	etaEast[simulationParameters.cells] = etaTemp[simulationParameters.cells] - hTemp[simulationParameters.cells] + hEastDownwind;
}

void uFaceValues(SimulationParameters simulationParameters, SolverParameters solverParameters, real* qFace, real* hFace, real* uFace)
{
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (hFace[i] <= solverParameters.tolDry)
		{
			uFace[i] = 0;
		}
		else
		{
			uFace[i] = qFace[i] / hFace[i];
		}
	}
}

void zStarIntermediateValues(SimulationParameters simulationParameters, real* etaWest, real* etaEast, real* hWest, real* hEast, real* zStarIntermediate)
{
	real a, b;
	
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		a = etaWest[i] - hWest[i];
		b = etaEast[i] - hEast[i];

		zStarIntermediate[i] = max(a, b);
	}
}

void deltaValues(SimulationParameters simulationParameters, real* etaFace, real* zStarIntermediate, real* deltaFace)
{
	real a;

	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		a = etaFace[i] - zStarIntermediate[i];
		deltaFace[i] = max(C(0.0), -a);
	}
}

// identical for both east and weat
void hStarValues(SimulationParameters simulationParameters, real* etaFaceValue, real* zStarIntermediate, real* hStar)
{
	real a;

	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		a = etaFaceValue[i] - zStarIntermediate[i];
		hStar[i] = max(C(0.0), a);
	}
}

// identical for both east and west
void qStarValues(SimulationParameters simulationParameters, real* uFaceValue, real* hStar, real* qStar)
{
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		qStar[i] = uFaceValue[i] * hStar[i];
	}
}

void zStarValues(SimulationParameters simulationParameters, real* deltaWest, real* deltaEast, real* zStarIntermediate, real* zStar)
{
	// wherever east =/= 0, west must be 0 and vice versa, so this works
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		zStar[i] = zStarIntermediate[i] - deltaEast[i] - deltaWest[i];
	}
}

void fluxHLL(SimulationParameters simulationParameters, SolverParameters solverParameters, real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* uWest, real* uEast, real* massFlux, real* momentumFlux)
{
	real aL, aR, hStar, uStar, aStar, sL, sR, massFL, massFR, momentumFL, momentumFR;

	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (hWestStar[i] <= solverParameters.tolDry && hEastStar[i] <= solverParameters.tolDry)
		{
			massFlux[i] = 0;
			momentumFlux[i] = 0;
			continue;
		}

		aL = sqrt(solverParameters.g * hWestStar[i]);
		aR = sqrt(solverParameters.g * hEastStar[i]);

		hStar = pow(((aL + aR) / 2 + (uWest[i] - uEast[i]) / 4), 2) / solverParameters.g;

		uStar = (uWest[i] + uEast[i]) / 2 + aL - aR;

		aStar = sqrt(solverParameters.g * hStar);

		if (hWestStar[i] <= solverParameters.tolDry)
		{
			sL = uEast[i] - 2 * aR;
		}
		else
		{
			sL = min(uWest[i] - aL, uStar - aStar);
		}

		if (hEastStar[i] <= solverParameters.tolDry)
		{
			sR = uWest[i] + 2 * aL;
		}
		else
		{
			sR = max(uEast[i] + aR, uStar + aStar);
		}

		massFL = qWestStar[i];
		massFR = qEastStar[i];

		momentumFL = uWest[i] * qWestStar[i] + solverParameters.g / 2 * pow(hWestStar[i], 2);
		momentumFR = uEast[i] * qEastStar[i] + solverParameters.g / 2 * pow(hEastStar[i], 2);

		if (sL >= 0)
		{
			massFlux[i] = massFL;
			momentumFlux[i] = momentumFL;
		}
		else if (sL < 0 && sR >= 0)
		{
			massFlux[i] = (sR * massFL - sL * massFR + sL * sR * (hEastStar[i] - hWestStar[i])) / (sR - sL);
			momentumFlux[i] = (sR * momentumFL - sL * momentumFR + sL * sR * (qEastStar[i] - qWestStar[i])) / (sR - sL);
		}
		else if (sR < 0)
		{
			massFlux[i] = massFR;
			momentumFlux[i] = momentumFR;
		}
	}
}

void hBarValues(SimulationParameters simulationParameters, real* hWestStar, real* hEastStar, real* hBar)
{
	// essentially 0th order projection but taking into account east/west locality
	for (int i = 0; i < simulationParameters.cells; i++)
	{
		hBar[i] = (hEastStar[i] + hWestStar[i + 1]) / 2;
	}
}

void momentumFV1OperatorValues(SimulationParameters simulationParameters, SolverParameters solverParameters, real dx, real* hBar, real* zBar, real* momentumFlux, real* momentumFV1Operator)
{
	real a = -(1 / dx);
	real b, c;

	for (int i = 0; i < simulationParameters.cells; i++)
	{
		b = momentumFlux[i + 1] - momentumFlux[i];
		c = 2 * sqrt(C(3.0)) * solverParameters.g * hBar[i] * zBar[i];

		momentumFV1Operator[i] = a * (b + c);
	}
}

void frictionImplicit(SimulationParameters simulationParameters, SolverParameters solverParameters, real dt, real* hWithBC, real* qWithBC)
{
	real qf, hf, u, Sf, D, Cf;
	
	for (int i = 0; i < simulationParameters.cells + 2; i++)
	{
		qf = qWithBC[i];
		hf = hWithBC[i];

		if (hWithBC[i] > solverParameters.tolDry && abs(qf) > solverParameters.tolDry)
		{
			u = qf / hf;

			Cf = solverParameters.g * pow(simulationParameters.manning, C(2.0)) / pow(hf, C(1.0)/C(3.0));

			Sf = -Cf * abs(u) * u;

			D = 1 + 2 * dt * Cf * abs(u) / hf;

			// Update
			qWithBC[i] += dt * Sf / D;
		}
	}
}