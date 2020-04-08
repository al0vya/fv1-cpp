#include "uFaceValues.h"

real uFaceValues(SolverParameters solverParameters, real qFace, real hFace)
{
	if (hFace <= solverParameters.tolDry)
	{
		return 0;
	}
	else
	{
		return qFace / hFace;
	}
}