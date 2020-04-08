#pragma once

#include <cmath>

#include "real.h"
#include "SolverParameters.h"
#include "SimulationParameters.h"

real frictionImplicit(SimulationParameters simulationParameters, SolverParameters solverParameters, real dt, real h, real q);