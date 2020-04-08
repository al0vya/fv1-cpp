#include "fluxHLL.h"

void fluxHLL(SimulationParameters simulationParameters, SolverParameters solverParameters, real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* uWest, real* uEast, real* massFlux, real* momentumFlux)
{
	real aL, aR, hStar, uStar, aStar, sL, sR, massFL, massFR, momentumFL, momentumFR;

	for (int i = 0; i < simulationParameters.cells + 1; i++)
	{
		if (hWestStar[i] <= solverParameters.tolDry && hEastStar[i] <= solverParameters.tolDry)
		{
			massFlux[i] = 0;
			momentumFlux[i] = 0;
			continue;
		}

		aL = sqrt(solverParameters.g * hWestStar[i]);
		aR = sqrt(solverParameters.g * hEastStar[i]);

		hStar = pow(((aL + aR) / 2 + (uWest[i] - uEast[i]) / 4), 2) / solverParameters.g;

		uStar = (uWest[i] + uEast[i]) / 2 + aL - aR;

		aStar = sqrt(solverParameters.g * hStar);

		if (hWestStar[i] <= solverParameters.tolDry)
		{
			sL = uEast[i] - 2 * aR;
		}
		else
		{
			sL = min(uWest[i] - aL, uStar - aStar);
		}

		if (hEastStar[i] <= solverParameters.tolDry)
		{
			sR = uWest[i] + 2 * aL;
		}
		else
		{
			sR = max(uEast[i] + aR, uStar + aStar);
		}

		massFL = qWestStar[i];
		massFR = qEastStar[i];

		momentumFL = uWest[i] * qWestStar[i] + solverParameters.g / 2 * pow(hWestStar[i], 2);
		momentumFR = uEast[i] * qEastStar[i] + solverParameters.g / 2 * pow(hEastStar[i], 2);

		if (sL >= 0)
		{
			massFlux[i] = massFL;
			momentumFlux[i] = momentumFL;
		}
		else if (sL < 0 && sR >= 0)
		{
			massFlux[i] = (sR * massFL - sL * massFR + sL * sR * (hEastStar[i] - hWestStar[i])) / (sR - sL);
			momentumFlux[i] = (sR * momentumFL - sL * momentumFR + sL * sR * (qEastStar[i] - qWestStar[i])) / (sR - sL);
		}
		else if (sR < 0)
		{
			massFlux[i] = massFR;
			momentumFlux[i] = momentumFR;
		}
	}
}