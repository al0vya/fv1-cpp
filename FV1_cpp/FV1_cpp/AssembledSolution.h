#pragma once

#include "real.h"

typedef struct AssembledSolution
{
	real* q_BC;
	real* h_BC;
	real* z_BC;
	int length;
	bool is_copy = false;

	AssembledSolution(const int& num_cells)
	:
		length(num_cells)
	{
		q_BC = new real[num_cells + 2];
		h_BC = new real[num_cells + 2];
		z_BC = new real[num_cells + 2];
	}

	AssembledSolution(const AssembledSolution& original) { is_copy = true; *this = original; }

	~AssembledSolution()
	{
		if (!is_copy)
		{
			delete[] q_BC;
			delete[] h_BC;
			delete[] z_BC;
		}
	}

} AssembledSolution;