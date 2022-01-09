#pragma once

#include "real.h"

typedef struct NodalValues
{
	real* q;
	real* h;
	real* z;
	real* x;
	bool is_copy = false;

	NodalValues(const int& num_interfaces)
	{
		q = new real[num_interfaces];
		h = new real[num_interfaces];
		z = new real[num_interfaces];
		x = new real[num_interfaces];
	}

	NodalValues(const NodalValues& original) { is_copy = true; *this = original; }

	~NodalValues()
	{
		if (!is_copy)
		{
			delete[] q;
			delete[] h;
			delete[] z;
			delete[] x;
		}
	}

} NodalValues;