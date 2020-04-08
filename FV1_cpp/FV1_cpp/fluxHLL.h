#pragma once
#include <cmath>
#include <algorithm>

using namespace std;

#include "real.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"

void fluxHLL(SimulationParameters simulationParameters, SolverParameters solverParameters, real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* uWest, real* uEast, real* massFlux, real* momentumFlux);