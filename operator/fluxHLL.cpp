#include "fluxHLL.h"

void fluxHLL
(
	SimulationParameters& simulationParameters, 
	SolverParameters&     solverParameters,
	StarValues&           star_vals,
	Fluxes&               fluxes
)
{
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (star_vals.h_west[i] <= solverParameters.tol_dry && star_vals.h_east[i] <= solverParameters.tol_dry)
		{
			fluxes.mass[i] = 0;
			fluxes.momentum[i] = 0;
		}
		else
		{
			real uWest = (star_vals.h_west[i] <= solverParameters.tol_dry) ? 0 : star_vals.q_west[i] / star_vals.h_west[i];
			real uEast = (star_vals.h_east[i] <= solverParameters.tol_dry) ? 0 : star_vals.q_east[i] / star_vals.h_east[i];

			real aWest = std::sqrt(solverParameters.g * star_vals.h_west[i]);
			real aEast = std::sqrt(solverParameters.g * star_vals.h_east[i]);

			real hStar = pow(((aWest + aEast) / 2 + (uWest - uEast) / 4), 2) / solverParameters.g;

			real uStar = (uWest + uEast) / 2 + aWest - aEast;

			real aStar = std::sqrt(solverParameters.g * hStar);

			real sWest = (star_vals.h_west[i] <= solverParameters.tol_dry) ? uEast - 2 * aEast : std::fmin(uWest - aWest, uStar - aStar);
			real sEast = (star_vals.h_east[i] <= solverParameters.tol_dry) ? uWest + 2 * aWest : std::fmax(uEast + aEast, uStar + aStar);

			real massWest = star_vals.q_west[i];
			real massEast = star_vals.q_east[i];

			real momentumWest = uWest * star_vals.q_west[i] + 0.5 * solverParameters.g * pow(star_vals.h_west[i], 2);
			real momentumEast = uEast * star_vals.q_east[i] + 0.5 * solverParameters.g * pow(star_vals.h_east[i], 2);

			if (sWest >= 0)
			{
				fluxes.mass[i] = massWest;
				fluxes.momentum[i] = momentumWest;
			}
			else if (sWest < 0 && sEast >= 0)
			{
				fluxes.mass[i] = (sEast * massWest - sWest * massEast + sWest * sEast * (star_vals.h_east[i] - star_vals.h_west[i])) / (sEast - sWest);
				fluxes.momentum[i] = (sEast * momentumWest - sWest * momentumEast + sWest * sEast * (star_vals.q_east[i] - star_vals.q_west[i])) / (sEast - sWest);
			}
			else if (sEast < 0)
			{
				fluxes.mass[i] = massEast;
				fluxes.momentum[i] = momentumEast;
			}
		}
	}

}