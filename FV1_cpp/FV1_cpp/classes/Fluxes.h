#pragma once

#include "../classes/real.h"

typedef struct Fluxes
{
	real* mass;
	real* momentum;

	bool is_copy = false;

	Fluxes(const int& num_interfaces)
	{
		mass     = new real[num_interfaces];
		momentum = new real[num_interfaces];
	}

	Fluxes(const Fluxes& original) { *this = original; is_copy = true; }

	~Fluxes()
	{
		if (!is_copy)
		{
			delete[] mass;
			delete[] momentum;
		}
	}

} Fluxes;