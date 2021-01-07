#pragma once

#include "real.h"

typedef struct
{
	real* qWithBC;
	real* hWithBC;
	real* zWithBC;
	real* dxLocalWithBC;
	real* x;
	int* activeIndices;
	int length;

} AssembledSolution;