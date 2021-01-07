#pragma once

#include "real.h"

typedef struct
{
	real* q_east;
	real* h_east;
	real* z_east;

	real* q_west;
	real* h_west;
	real* z_west;

} StarValues;