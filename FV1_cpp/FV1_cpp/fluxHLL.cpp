#include "fluxHLL.h"

void fluxHLL(SimulationParameters simulationParameters, SolverParameters solverParameters, real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* massFlux, real* momentumFlux)
{
	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (hWestStar[i] <= solverParameters.tolDry && hEastStar[i] <= solverParameters.tolDry)
		{
			massFlux[i] = 0;
			momentumFlux[i] = 0;
		}
		else
		{
			real uWest = (hWestStar[i] <= solverParameters.tolDry) ? 0 : qWestStar[i] / hWestStar[i];
			real uEast = (hEastStar[i] <= solverParameters.tolDry) ? 0 : qEastStar[i] / hEastStar[i];

			real aWest = sqrt(solverParameters.g * hWestStar[i]);
			real aEast = sqrt(solverParameters.g * hEastStar[i]);

			real hStar = pow(((aWest + aEast) / 2 + (uWest - uEast) / 4), 2) / solverParameters.g;

			real uStar = (uWest + uEast) / 2 + aWest - aEast;

			real aStar = sqrt(solverParameters.g * hStar);

			real sWest = (hWestStar[i] <= solverParameters.tolDry) ? uEast - 2 * aEast : min(uWest - aWest, uStar - aStar);
			real sEast = (hEastStar[i] <= solverParameters.tolDry) ? uWest + 2 * aWest : max(uEast + aEast, uStar - aStar);

			real massWest = qWestStar[i];
			real massEast = qEastStar[i];

			real momentumWest = uWest * qWestStar[i] + 0.5 * solverParameters.g * pow(hWestStar[i], 2);
			real momentumEast = uEast * qEastStar[i] + 0.5 * solverParameters.g * pow(hEastStar[i], 2);

			if (sWest >= 0)
			{
				massFlux[i] = massWest;
				momentumFlux[i] = momentumWest;
			}
			else if (sWest < 0 && sEast >= 0)
			{
				massFlux[i] = (sEast * massWest - sWest * massEast + sWest * sEast * (hEastStar[i] - hWestStar[i])) / (sEast - sWest);
				momentumFlux[i] = (sEast * momentumWest - sWest * momentumEast + sWest * sEast * (qEastStar[i] - qWestStar[i])) / (sEast - sWest);
			}
			else if (sEast < 0)
			{
				massFlux[i] = massEast;
				momentumFlux[i] = momentumEast;
			}
		}
	}

}