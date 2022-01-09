#pragma once

#include <cmath>

#include "real.h"
#include "SolverParameters.h"
#include "SimulationParameters.h"

real friction_implicit(SimulationParameters simulationParameters, SolverParameters solverParameters, real dt, real h, real q);