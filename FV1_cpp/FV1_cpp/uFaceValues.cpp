#include "uFaceValues.h"

real uFaceValues(SolverParameters solverParameters, real qFace, real hFace)
{
	if (hFace <= solverParameters.tol_dry)
	{
		return 0;
	}
	else
	{
		return qFace / hFace;
	}
}