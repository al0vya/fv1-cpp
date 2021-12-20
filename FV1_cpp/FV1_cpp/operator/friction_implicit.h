#pragma once

#include <cmath>

#include "../classes/real.h"
#include "../classes/SolverParameters.h"
#include "../classes/SimulationParameters.h"

real friction_implicit(SimulationParameters simulationParameters, SolverParameters solverParameters, real dt, real h, real q);