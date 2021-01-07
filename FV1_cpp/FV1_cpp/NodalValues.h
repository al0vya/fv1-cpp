#pragma once

#include "real.h"

typedef struct NodalValues
{
	real* q;
	real* h;
	real* z;
	real* x;

} NodalValues;