#include "real.h"
#include "baselineMesh.h"

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