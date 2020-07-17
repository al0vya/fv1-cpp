#include "fluxHLLNew.h"

void fluxHLLNew(SimulationParameters simulationParameters, SolverParameters solverParameters, real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* massFlux, real* momentumFlux)
{
	real aL, aR, hStar, uStar, aStar, sL, sR, massFL, massFR, momentumFL, momentumFR, uWest, uEast;

	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (hWestStar[i] <= solverParameters.tolDry && hEastStar[i] <= solverParameters.tolDry)
		{
			massFlux[i] = 0;
			momentumFlux[i] = 0;
			continue;
		}

		uWest = (hWestStar[i] <= solverParameters.tolDry) ? 0 : qWestStar[i] / hWestStar[i];
		uEast = (hEastStar[i] <= solverParameters.tolDry) ? 0 : qEastStar[i] / hEastStar[i];

		aL = sqrt(solverParameters.g * hWestStar[i]);
		aR = sqrt(solverParameters.g * hEastStar[i]);

		hStar = pow(((aL + aR) / 2 + (uWest - uEast) / 4), 2) / solverParameters.g;

		uStar = (uWest + uEast) / 2 + aL - aR;

		aStar = sqrt(solverParameters.g * hStar);

		if (hWestStar[i] <= solverParameters.tolDry)
		{
			sL = uEast - 2 * aR;
		}
		else
		{
			sL = min(uWest - aL, uStar - aStar);
		}

		if (hEastStar[i] <= solverParameters.tolDry)
		{
			sR = uWest + 2 * aL;
		}
		else
		{
			sR = max(uEast + aR, uStar + aStar);
		}

		massFL = qWestStar[i];
		massFR = qEastStar[i];

		momentumFL = uWest * qWestStar[i] + solverParameters.g / 2 * pow(hWestStar[i], 2);
		momentumFR = uEast * qEastStar[i] + solverParameters.g / 2 * pow(hEastStar[i], 2);

		if (sL >= 0)
		{
			massFlux[i] = massFL;
			momentumFlux[i] = momentumFL;
		}
		else if (sL < 0 && sR >= 0)
		{
			massFlux[i] = (sR * massFL - sL * massFR + sL * sR * (hEastStar[i] - hWestStar[i])) / (sR - sL);
			momentumFlux[i] = (sR * momentumFL - sL * momentumFR + sL * sR * (qEastStar[i] - qWestStar[i])) / (sR - sL);
			int dummy = 1;
		}
		else if (sR < 0)
		{
			massFlux[i] = massFR;
			momentumFlux[i] = momentumFR;
		}
	}
}