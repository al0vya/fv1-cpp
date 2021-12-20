#pragma once
#include "../classes/real.h"

typedef struct SolverParameters
{
	real CFL;
	real tol_dry;
	real g;

} SolverParameters;