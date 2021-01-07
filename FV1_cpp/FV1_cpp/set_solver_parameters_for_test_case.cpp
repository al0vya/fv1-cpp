#include "set_solver_parameters_for_test_case.h"
#include "real.h"

SolverParameters set_solver_parameters_for_test_case()
{
	SolverParameters solverParameters;

	solverParameters.CFL = C(0.33);
	solverParameters.tol_dry = C(1e-4);
	solverParameters.g = C(9.80665);

	return solverParameters;
}