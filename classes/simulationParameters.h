#pragma once

#include "real.h"

typedef struct SimulationParameters
{
	int cells;

	real xmin;
	real xmax;

	real simulationTime;

	real manning;
} SimulationParameters;