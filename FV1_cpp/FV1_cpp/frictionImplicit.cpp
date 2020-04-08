#include "frictionImplicit.h"

real frictionImplicit(SimulationParameters simulationParameters, SolverParameters solverParameters, real dt, real h, real q)
{
	real u, Sf, D, Cf;

	if (h > solverParameters.tolDry && abs(q) > solverParameters.tolDry)
	{
		u = q / h;

		Cf = solverParameters.g * pow(simulationParameters.manning, C(2.0)) / pow(h, C(1.0) / C(3.0));

		Sf = -Cf * abs(u) * u;

		D = 1 + 2 * dt * Cf * abs(u) / h;

		// Update
		return dt * Sf / D;
	}
	else
	{
		return 0;
	}
}