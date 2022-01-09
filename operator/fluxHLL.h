#pragma once
#include <cmath>
#include <algorithm>

using namespace std;

#include "real.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "StarValues.h"
#include "Fluxes.h"

void fluxHLL
(
	SimulationParameters& simulationParameters, 
	SolverParameters&     solverParameters,
	StarValues&           star_vals,
	Fluxes&               fluxes
);