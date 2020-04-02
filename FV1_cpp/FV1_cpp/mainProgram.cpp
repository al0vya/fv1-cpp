#include <stdlib.h>
#include <stdio.h>


typedef float real;

// defining simulation parameters
int cells = 50;
int xmin = 0;
int xmax = 50;

real simulationTime = 2.5;

real g = (real)9.80665;
real manning = 0.0;

int hl = 6;
int hr = 2;

int ql = 0;
int qr = 0;

int reflectUp = 1;
int reflectDown = 1;

int hImposedUp = 0;
int qxImposedUp = 0;

int hImposedDown = 0;
int qxImposedDown = 0;

real timeNow = 0;

real CFL = (real)0.33;

real tolDry = (real)1e-3;

real dt = (real)1e-4;

// declaring helper functions
void baselineMesh(real dx, real* x, real* x_int);
void bedDataDamBreak(real* z_int);
void qInitialDamBreak(real* x_int, real* q_int);
void hInitialDamBreak(real* z_int, real* x_int, real* h_int);
void modalProjection(real* u_int, real* u);

int main()
{

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

	// test they've been initialised
	for (i = 0; i < cells+1; i++) {
		printf("%f, %f\n", x_int[i], h_int[i]);
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

// note this function works on the z NODES
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