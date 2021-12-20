#pragma once

#include "../classes/real.h"

typedef struct BarValues
{
	real* h;
	real* z;

	bool is_copy = false;

	BarValues(const int& num_cells)
	{
		h = new real[num_cells];
		z = new real[num_cells];
	}

	BarValues(const BarValues& original) { *this = original; is_copy = true; }

	~BarValues()
	{
		if (!is_copy)
		{
			delete[] h;
			delete[] z;
		}
	}

} BarValues;