#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
using namespace std;

typedef double real;

// Defining simulation parameters //
int cells = 10;
real xmin = 0;
real xmax = 50;

real simulationTime = 2.5;

real g = (real)9.80665;
real manning = 0.0;

real hl = 6;
real hr = 2;

real ql = 0;
real qr = 0;

real reflectUp = 1;
real reflectDown = 1;

real hImposedUp = 0;
real qxImposedUp = 0;

real hImposedDown = 0;
real qxImposedDown = 0;

real timeNow = 0;

real CFL = 0.33f;

real tolDry = 1e-3f;

real dt = 1e-4;

// declaring helper functions
void baselineMesh(real dx, real* x, real* x_int);
void bedDataDamBreak(real* z_int);
void qInitialDamBreak(real* x_int, real* q_int);
void hInitialDamBreak(real* z_int, real* x_int, real* h_int);
void modalProjectionZeroOrder(real* u_int, real* u);
void modalProjectionFirstOrder(real* u_int, real* u);
void qAddGhostBoundaryConditions(real* q, real* qWithBC);
void hAddGhostBoundaryConditions(real* h, real* hWithBC);
void zAddGhostBoundaryConditions(real* z, real* zWithBC);
void hWestFaceValues(real hWestUpwind, real* hTemp, real* hWest);
void qWestFaceValues(real qWestUpwind, real* qTemp, real* qWest);
void etaWestFaceValues(real hWestUpwind, real* hTemp, real* etaTemp, real* etaWest);
void hEastFaceValues(real hEastDownwind, real* hTemp, real* hEast);
void qEastFaceValues(real qEastDownwind, real* qTemp, real* qEast);
void etaEastFaceValues(real hEastDownwind, real* hTemp, real* etaTemp, real* etaEast);
void uWestFaceValues(real* qWest, real* hWest, real* uWest);
void uEastFaceValues(real* qEast, real* hEast, real* uEast);
void zStarIntermediateValues(real* etaWest, real* etaEast, real* hWest, real* hEast, real* zStarIntermediate);
void deltaWestValues(real* etaWest, real* zStarIntermediate, real* deltaWest);
void deltaEastValues(real* etaEast, real* zStarIntermediate, real* deltaEast);
void hStarValues(real* etaFaceValue, real* zStarIntermediate, real* hStar);
void qStarValues(real* uFaceValue, real* hStar, real* qStar);
void zStarValues(real* deltaWest, real* deltaEast, real* zStarIntermediate, real* zStar);
void fluxHLL(real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* uWest, real* uEast, real* massFlux, real* momentumFlux);
void hBarValues(real* hWestStar, real* hEastStar, real* hBar);
void massFV1OperatorValues(real dx, real* massFlux, real* massFV1Operator);
void momentumFV1OperatorValues(real dx, real* hBar, real* zBar, real* momentumFlux, real* momentumFV1Operator);

int main()
{
	// quintessential for-loop index
	int i;
	real u;
	real dtCFL;

	// coarsest cell size
	real dx = (xmax - xmin) / (real)cells;

	// allocate buffers for the cell centres and interfaces
	real* x = new real[cells];
	real* x_int = new real[cells + 1];

	// initialise baseline mesh
	baselineMesh(dx, x, x_int);

	// allocate buffers for flow modes q, h, z
	real* q = new real[cells];
	real* h = new real[cells];
	real* z = new real[cells];

	// allocate buffers for the corresponding interface Valuess
	real* q_int = new real[cells + 1];
	real* h_int = new real[cells + 1];
	real* z_int = new real[cells + 1];

	// initialise the interface Valuess
	bedDataDamBreak(z_int);
	qInitialDamBreak(x_int, q_int);
	hInitialDamBreak(z_int, x_int, h_int);

	// project to find the modes
	modalProjectionZeroOrder(q_int, q);
	modalProjectionZeroOrder(h_int, h);
	modalProjectionZeroOrder(z_int, z);

	// allocate buffers for modes with BCs
	real* qWithBC = new real[cells + 2];
	real* hWithBC = new real[cells + 2];
	real* zWithBC = new real[cells + 2];

	// allocate buffers for RK2 step
	real* qWithBCNew = new real[cells + 2];
	real* hWithBCNew = new real[cells + 2];

	real* qTemp = new real[cells + 2];
	real* hTemp = new real[cells + 2];
	real* zTemp = new real[cells + 2];

	// placeholder variables for substitution
	real a, b, c;

	// allocate true/false buffer for dry cells
	bool* dry = new bool[cells + 2];

	// placeholder variables for checking
	real hLocal, hBackward, hForward, hMax;

	real* etaTemp = new real[cells + 2];

	// allocating buffers for eastern and western interface Valuess
	real* qEast = new real[cells + 1];
	real* hEast = new real[cells + 1];
	real* etaEast = new real[cells + 1];

	real* qWest = new real[cells + 1];
	real* hWest = new real[cells + 1];
	real* etaWest = new real[cells + 1];

	// allocating buffers for positivity preserving nodes
	real* qEastStar = new real[cells + 1];
	real* hEastStar = new real[cells + 1];

	real* qWestStar = new real[cells + 1];
	real* hWestStar = new real[cells + 1];

	real* zStarIntermediate = new real[cells + 1];
	real* zStar = new real[cells + 1];

	real* uWest = new real[cells + 1];
	real* uEast = new real[cells + 1];

	real* deltaWest = new real[cells + 1];
	real* deltaEast = new real[cells + 1];

	// allocating buffers for numerical fluxes with HLL solver
	real* massFlux = new real[cells + 1];
	real* momentumFlux = new real[cells + 1];

	// allocating buffers for positivity preserving MODES
	real* hBar = new real[cells];
	real* zBar = new real[cells];

	// allocating buffer for fv1Operator values
	real* massFV1Operator = new real[cells];
	real* momentumFV1Operator = new real[cells];

	// WHILE LOOP STARTS FROM HERE //

	while (timeNow < simulationTime)
	{
		if (timeNow - simulationTime > 0)
		{
			timeNow -= dt;
			dt = simulationTime - timeNow;
			timeNow += dt;
		}

		// initialise modes with BCs
		qAddGhostBoundaryConditions(q, qWithBC);
		hAddGhostBoundaryConditions(h, hWithBC);
		zAddGhostBoundaryConditions(z, zWithBC);

		// extract upwind and downwind modes for use in fv1Operator
		real hWestUpwind = hWithBC[0];
		real qWestUpwind = qWithBC[0];

		real hEastDownwind = hWithBC[cells + 1];
		real qEastDownwind = qWithBC[cells + 1];

		// initialise RK2 buffers
		for (i = 0; i < cells + 2; i++)
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

		// ghost cells are always 0 i.e. false/wet
		dry[0] = false;
		dry[cells + 1] = false;

		// initialising dry vs wet cells, ignore ghost cells
		for (i = 1; i < cells + 1; i++)
		{
			hLocal = hTemp[i];
			hBackward = hTemp[i - 1];
			hForward = hTemp[i + 1];

			hMax = max(hLocal, hBackward);
			hMax = max(hForward, hMax);

			// dry[] hasn't been initialised so else statement also needed
			if (hMax <= tolDry)
			{
				dry[i] = true;
			}
			else
			{
				dry[i] = false;
			}
		}

		for (i = 0; i < cells + 1; i++)
		{
			etaTemp[i] = hTemp[i] + zTemp[i];
		}

		// initialising interface Valuess
		qEastFaceValues(qEastDownwind, qTemp, qEast);
		hEastFaceValues(hEastDownwind, hTemp, hEast);
		etaEastFaceValues(hEastDownwind, hTemp, etaTemp, etaEast);

		qWestFaceValues(qWestUpwind, qTemp, qWest);
		hWestFaceValues(hWestUpwind, hTemp, hWest);
		etaWestFaceValues(hWestUpwind, hTemp, etaTemp, etaWest);

		// initialising velocity interface values
		uWestFaceValues(qWest, hWest, uWest);
		uEastFaceValues(qEast, hEast, uEast);

		zStarIntermediateValues(etaWest, etaEast, hWest, hEast, zStarIntermediate);

		deltaWestValues(etaWest, zStarIntermediate, deltaWest);
		deltaEastValues(etaEast, zStarIntermediate, deltaEast);

		// initialising positivity preserving nodes
		hStarValues(etaWest, zStarIntermediate, hWestStar);
		hStarValues(etaEast, zStarIntermediate, hEastStar);

		qStarValues(uWest, hWestStar, qWestStar);
		qStarValues(uEast, hEastStar, qEastStar);

		zStarValues(deltaWest, deltaEast, zStarIntermediate, zStar);

		// initialising numerical fluxes
		fluxHLL(hWestStar, hEastStar, qWestStar, qEastStar, uWest, uEast, massFlux, momentumFlux);

		hBarValues(hWestStar, hEastStar, hBar);

		modalProjectionFirstOrder(zStar, zBar);

		massFV1OperatorValues(dx, massFlux, massFV1Operator);
		momentumFV1OperatorValues(dx, hBar, zBar, momentumFlux, momentumFV1Operator);

		// FV1 operator increment, skip ghosts cells
		for (i = 1; i < cells + 1; i++)
		{
			// skip increment in dry cells
			if (dry[i])
			{
				hWithBCNew[i] = hWithBC[i];
				qWithBCNew[i] = qWithBC[i];
			}
			else
			{
				hWithBCNew[i] = hWithBC[i] + dt * massFV1Operator[i - 1];
				qWithBCNew[i] = qWithBC[i] + dt * momentumFV1Operator[i - 1];
			}

			if (hWithBCNew[i] <= tolDry)
			{
				qWithBCNew[i] = 0;
			}
		}

		// update and CFL time step adjustment
		dt = 1e9;

		for (i = 0; i < cells; i++)
		{
			h[i] = hWithBCNew[i + 1];
			q[i] = qWithBCNew[i + 1];

			if (h[i] <= tolDry)
			{
				continue;
			}
			else
			{
				u = q[i] / h[i];
				dtCFL = CFL * dx / (abs(u) + sqrt(g * h[i]));
				dt = min(dt, dtCFL);
			}
		}

		printf("Updated time step, h value: %f, q value %f.\n", h[0], q[0]);
	}

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

void baselineMesh(real dx, real* x, real* x_int)
{
	int i;
	real a, b;

	x_int[0] = xmin;
	x_int[cells] = xmax;

	for (i = 1; i < cells; i++) {
		a = x_int[i - 1];
		b = a + dx;

		x_int[i] = b;
		x[i - 1] = (a + b) / 2;
	}

	// final cell centre since it's not done in the for loop
	x[cells - 1] = (x_int[cells] + x_int[cells - 1]) / 2;
}

// note this function works on the z NODES i.e. interfaces
void bedDataDamBreak(real* z_int)
{
	int i;

	// for (cells+1) interfaces
	for (i = 0; i < cells + 1; i++)
	{
		z_int[i] = 0;
	}
}

// note this function works on the interfaces
void qInitialDamBreak(real* x_int, real* q_int)
{
	int i;

	// recall for loop is for (cells+1) interfaces
	for (i = 0; i < cells + 1; i++)
	{
		if (x_int[i] <= 32.5)
		{
			q_int[i] = ql;
		}
		else
		{
			q_int[i] = qr;
		}
	}
}

void hInitialDamBreak(real* z_int, real* x_int, real* h_int)
{
	real etaLeft = hl;
	real etaRight = hr;

	int i;

	// recall for loop is for (cells+1) interfaces
	for (i = 0; i < cells + 1; i++)
	{
		if (x_int[i] <= 25)
		{
			h_int[i] = etaLeft - z_int[i];
		}
		else
		{
			h_int[i] = etaRight - z_int[i];
		}
	}
}

void modalProjectionZeroOrder(real* u_int, real* u)
{
	int i;

	// for modes, so loops overs cells rather than (cells+1) interfaces
	for (i = 0; i < cells; i++)
	{
		u[i] = (u_int[i] + u_int[i + 1]) / 2;
	}
}

void modalProjectionFirstOrder(real* u_int, real* u)
{
	for (int i = 0; i < cells; i++)
	{
		u[i] = (u_int[i + 1] - u_int[i]) / 2 * sqrt(3);
	}
}

void qAddGhostBoundaryConditions(real* q, real* qWithBC)
{
	real qUp = reflectUp * q[0];

	if (qxImposedUp > 0)
	{
		qUp = qxImposedUp;
	}

	real qDown = reflectDown * q[cells - 1];

	if (qxImposedDown > 0)
	{
		qDown = qxImposedDown;
	}

	qWithBC[0] = qUp;

	// there are cells + 2 elements inc BCs
	qWithBC[cells + 1] = qDown;

	int i;

	for (i = 1; i < cells + 1; i++)
	{
		qWithBC[i] = q[i - 1];
	}
}

void hAddGhostBoundaryConditions(real* h, real* hWithBC)
{
	real hUp = h[0];

	if (hImposedUp > 0)
	{
		hUp = hImposedUp;
	}

	real hDown = h[cells - 1];

	if (hImposedDown > 0)
	{
		hDown = hImposedDown;
	}

	hWithBC[0] = hUp;

	// there are cells + 2 elements inc BCs
	hWithBC[cells + 1] = hDown;

	int i;

	for (i = 1; i < cells + 1; i++)
	{
		hWithBC[i] = h[i - 1];
	}
}

void zAddGhostBoundaryConditions(real* z, real* zWithBC)
{
	real zUp = z[0];
	real zDown = z[cells - 1];

	zWithBC[0] = zUp;

	// there are cells + 2 elements inc BCs
	zWithBC[cells + 1] = zDown;

	int i;

	for (i = 1; i < cells + 1; i++)
	{
		zWithBC[i] = z[i - 1];
	}
}

void qLocalFaceValues(int i, real* qTemp, real* q)
{
	q[i] = qTemp[i];
}

void hLocalFaceValues(int i, real* hTemp, real* h)
{
	h[i] = hTemp[i];
}

void etaLocalFaceValues(int i, real* hTemp, real* z, real* eta)
{
	eta[i] = hTemp[i] + z[i];
}

void hWestFaceValues(real hWestUpwind, real* hTemp, real* hWest)
{
	hWest[0] = hWestUpwind;

	int i;

	// start from 1 on hTemp to avoid ghost cell
	for (i = 1; i < cells + 1; i++)
	{
		hWest[i] = hTemp[i];
	}
}

void qWestFaceValues(real qWestUpwind, real* qTemp, real* qWest)
{
	qWest[0] = qWestUpwind;

	int i;

	// start from 1 on qTemp to avoid ghost cell
	for (i = 1; i < cells + 1; i++)
	{
		qWest[i] = qTemp[i];
	}
}

void etaWestFaceValues(real hWestUpwind, real* hTemp, real* etaTemp, real* etaWest)
{
	// indexing shifted by one, since uTemp[0] is the ghost mode
	etaWest[0] = etaTemp[1] - hTemp[1] + hWestUpwind;

	int i;

	// similarly shifted, start from 1 on temp to avoid ghost cell
	for (i = 1; i < cells + 1; i++)
	{
		etaWest[i] = etaTemp[i];
	}
}

void hEastFaceValues(real hEastDownwind, real* hTemp, real* hEast)
{
	int i;

	// [i+1] on hTemp to skip ghost cell
	for (i = 0; i < cells; i++)
	{
		hEast[i] = hTemp[i + 1];
	}

	hEast[cells] = hEastDownwind;
}

void qEastFaceValues(real qEastDownwind, real* qTemp, real* qEast)
{
	int i;

	// [i+1] on qTemp to skip ghost cell
	for (i = 0; i < cells; i++)
	{
		qEast[i] = qTemp[i + 1];
	}

	qEast[cells] = qEastDownwind;
}

void etaEastFaceValues(real hEastDownwind, real* hTemp, real* etaTemp, real* etaEast)
{
	int i;

	for (i = 0; i < cells; i++)
	{
		etaEast[i] = etaTemp[i + 1];
	}

	etaEast[cells] = etaTemp[cells] - hTemp[cells] + hEastDownwind;
}

void uWestFaceValues(real* qWest, real* hWest, real* uWest)
{
	for (int i = 0; i < cells + 1; i++)
	{
		if (hWest[i] <= tolDry)
		{
			uWest[i] = 0;
		}
		else
		{
			uWest[i] = qWest[i] / hWest[i];
		}
	}
}

void uEastFaceValues(real* qEast, real* hEast, real* uEast)
{
	for (int i = 0; i < cells + 1; i++)
	{
		if (hEast[i] <= tolDry)
		{
			uEast[i] = 0;
		}
		else
		{
			uEast[i] = qEast[i] / hEast[i];
		}
	}
}

void zStarIntermediateValues(real* etaWest, real* etaEast, real* hWest, real* hEast, real* zStarIntermediate)
{
	real a, b;
	
	for (int i = 0; i < cells + 1; i++)
	{
		a = etaWest[i] - hWest[i];
		b = etaEast[i] - hEast[i];

		zStarIntermediate[i] = max(a, b);
	}
}

void deltaWestValues(real* etaWest, real* zStarIntermediate, real* deltaWest)
{
	real a;

	for (int i = 0; i < cells + 1; i++)
	{
		a = etaWest[i] - zStarIntermediate[i];
		deltaWest[i] = max(0.0, -a);
	}
}

void deltaEastValues(real* etaEast, real* zStarIntermediate, real* deltaEast)
{
	real a;

	for (int i = 0; i < cells + 1; i++)
	{
		a = etaEast[i] - zStarIntermediate[i];
		deltaEast[i] = max(0.0, -a);
	}
}

// identical for both east and weat
void hStarValues(real* etaFaceValue, real* zStarIntermediate, real* hStar)
{
	real a;

	for (int i = 0; i < cells + 1; i++)
	{
		a = etaFaceValue[i] - zStarIntermediate[i];
		hStar[i] = max(0.0, a);
	}
}

// identical for both east and west
void qStarValues(real* uFaceValue, real* hStar, real* qStar)
{
	for (int i = 0; i < cells + 1; i++)
	{
		qStar[i] = uFaceValue[i] * hStar[i];
	}
}

void zStarValues(real* deltaWest, real* deltaEast, real* zStarIntermediate, real* zStar)
{
	// wherever east =/= 0, west must be 0 and vice versa, so this works
	for (int i = 0; i < cells + 1; i++)
	{
		zStar[i] = zStarIntermediate[i] - deltaEast[i] - deltaWest[i];
	}
}

void fluxHLL(real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* uWest, real* uEast, real* massFlux, real* momentumFlux)
{
	real aL, aR, hStar, uStar, aStar, sL, sR, massFL, massFR, momentumFL, momentumFR;

	for (int i = 0; i < cells + 1; i++)
	{
		if (hWestStar[i] <= tolDry && hEastStar[i] <= tolDry)
		{
			massFlux[i] = 0;
			momentumFlux[i] = 0;
			continue;
		}

		aL = sqrt(g * hWestStar[i]);
		aR = sqrt(g * hEastStar[i]);

		hStar = pow(((aL + aR) / 2 + (uWest[i] - uEast[i]) / 4), 2) / g;

		uStar = (uWest[i] + uEast[i]) / 2 + aL - aR;

		aStar = sqrt(g * hStar);

		if (hWestStar[i] <= tolDry)
		{
			sL = uEast[i] - 2 * aR;
		}
		else
		{
			sL = min(uWest[i] - aL, uStar - aStar);
		}

		if (hEastStar[i] <= tolDry)
		{
			sR = uWest[i] + 2 * aL;
		}
		else
		{
			sR = max(uEast[i] + aR, uStar + aStar);
		}

		massFL = qWestStar[i];
		massFR = qEastStar[i];

		momentumFL = uWest[i] * qWestStar[i] + g / 2 * pow(hWestStar[i], 2);
		momentumFR = uEast[i] * qEastStar[i] + g / 2 * pow(hEastStar[i], 2);

		if (sL >= 0)
		{
			massFlux[i] = massFL;
			momentumFlux[i] = momentumFL;
		}
		else if (sL < 0 && sR >= 0)
		{
			massFlux[i] = sR * massFL - sL * massFR + sL * sR * (hEastStar[i] - hWestStar[i]) / (sR - sL);
			momentumFlux[i] = sR * momentumFL - sL * momentumFR + sL * sR * (qEastStar[i] - qWestStar[i]) / (sR - sL);
		}
		else if (sR < 0)
		{
			massFlux[i] = massFR;
			momentumFlux[i] = momentumFR;
		}
	}
}

void hBarValues(real* hWestStar, real* hEastStar, real* hBar)
{
	// essentially 0th order projection but taking into account east/west locality
	for (int i = 0; i < cells; i++)
	{
		hBar[i] = (hEastStar[i] + hWestStar[i + 1]) / 2;
	}
}

void massFV1OperatorValues(real dx, real* massFlux, real* massFV1Operator)
{
	real a = -(1 / dx);

	for (int i = 0; i < cells; i++)
	{
		massFV1Operator[i] = a * (massFlux[i + 1] - massFlux[i]);
	}
}

void momentumFV1OperatorValues(real dx, real* hBar, real* zBar, real* momentumFlux, real* momentumFV1Operator)
{
	real a = -(1 / dx);

	for (int i = 0; i < cells; i++)
	{
		momentumFV1Operator[i] = a * (momentumFlux[i + 1] - momentumFlux[i]) + 2 * sqrt(3) * g * hBar[i] * zBar[i];
	}
}