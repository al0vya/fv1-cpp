#pragma once

#include "real.h"

typedef struct StarValues
{
	real* q_east;
	real* h_east;
	real* z_east;

	real* q_west;
	real* h_west;
	real* z_west;

	bool is_copy = false;

	StarValues(const int& num_interfaces)
	{
		q_east = new real[num_interfaces];
	    h_east = new real[num_interfaces];
	    z_east = new real[num_interfaces];
	    
	    q_west = new real[num_interfaces];
	    h_west = new real[num_interfaces];
	    z_west = new real[num_interfaces];
	}

	StarValues(const StarValues& original) { is_copy = true; *this = original; }

	~StarValues()
	{
		if (!is_copy)
		{
			delete[] q_east;
			delete[] h_east;
			delete[] z_east;

			delete[] q_west;
			delete[] h_west;
			delete[] z_west;
		}
	}

} StarValues;