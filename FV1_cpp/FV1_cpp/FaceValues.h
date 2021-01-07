#pragma once

#include "real.h"

typedef struct
{
	real* q_east;
	real* h_east;
	real* eta_east;

	real* q_west;
	real* h_west;
	real* eta_west;

} FaceValues;