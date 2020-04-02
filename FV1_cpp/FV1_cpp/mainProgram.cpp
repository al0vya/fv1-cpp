#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
using namespace std;

typedef double real;

// Defining simulation parameters //
int cells = 20;
real xmin = 0;
real xmax = 50;

real simulationTime = 2.5;

real g = (real)9.80665;
real manning = 0.0;

int hl = 6;
int hr = 2;

real ql = 0;
real qr = 0;

int reflectUp = 1;
int reflectDown = 1;

int hImposedUp = 0;
int qxImposedUp = 0;

int hImposedDown = 0;
int qxImposedDown = 0;

real timeNow = 0;

real CFL = 0.33f;

real tolDry = 1e-3f;

real dt = 1e-4f;

// declaring helper functions
void baselineMesh(real dx, real* x, real* x_int);
void bedDataDamBreak(real* z_int);
void qInitialDamBreak(real* x_int, real* q_int);
void hInitialDamBreak(real* z_int, real* x_int, real* h_int);
void modalProjection(real* u_int, real* u);
void qAddGhostBoundaryConditions(real* q, real* qWithBC);
void hAddGhostBoundaryConditions(real* h, real* hWithBC);
void zAddGhostBoundaryConditions(real* z, real* zWithBC);

int main()
{
	// quintessential for-loop index
	int i;

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

	// allocate buffers for the corresponding interfacial values
	real* q_int = new real[cells + 1];
	real* h_int = new real[cells + 1];
	real* z_int = new real[cells + 1];

	// initialise the interfacial values
	bedDataDamBreak(z_int);
	qInitialDamBreak(x_int, q_int);
	hInitialDamBreak(z_int, x_int, h_int);

	// project to find the modes
	modalProjection(q_int, q);
	modalProjection(h_int, h);
	modalProjection(z_int, z);

	// allocate buffers for modes with BCs
	real* qWithBC = new real[cells + 2];
	real* hWithBC = new real[cells + 2];
	real* zWithBC = new real[cells + 2];

	// initialise modes with BCs
	qAddGhostBoundaryConditions(q, qWithBC);
	hAddGhostBoundaryConditions(h, hWithBC);
	zAddGhostBoundaryConditions(z, zWithBC);

	// allocate buffers for RK2 step
	real* qWithBCNew = new real[cells + 2];
	real* hWithBCNew = new real[cells + 2];

	real* qTemp = new real[cells + 2];
	real* hTemp = new real[cells + 2];
	real* zTemp = new real[cells + 2];

	// placeholder variables for substitution
	real a, b, c;

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

	// allocate true/false buffer for dry cells
	bool* dry = new bool[cells + 2];

	// ghost cells are always 0 i.e. false
	dry[0] = false;
	dry[cells + 1] = false;

	// placeholder variables for checking
	real hLocal, hBackward, hForward, hMax;
	
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

	int j = 0;

	// test they've been initialised
	for (i = 0; i < cells+2; i++) {
		j++;
		printf("%f, %s\n", hTemp[i], dry[i] ? "true" : "false");
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
	int etaLeft = hl;
	int etaRight = hr;

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

void modalProjection(real* u_int, real* u)
{
	int i;

	// for modes, so loops overs cells rather than (cells+1) interfaces
	for (i = 0; i < cells; i++)
	{
		u[i] = (u_int[i] + u_int[i + 1]) / 2;
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
		hUp = hImposedDown;
	}

	real hDown = h[cells - 1];

	if (qxImposedDown > 0)
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

