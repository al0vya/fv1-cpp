#pragma once

#include "../classes/real.h"

typedef struct FaceValues
{
	real* q_east;
	real* h_east;
	real* eta_east;

	real* q_west;
	real* h_west;
	real* eta_west;

	bool is_copy = false;

	FaceValues(const int num_interfaces)
	{
		q_east   = new real[num_interfaces];
	    h_east   = new real[num_interfaces];
	    eta_east = new real[num_interfaces];
	    
	    q_west   = new real[num_interfaces];
	    h_west   = new real[num_interfaces];
	    eta_west = new real[num_interfaces];
	}

	FaceValues(const FaceValues& original) { is_copy = true; *this = original; }

	~FaceValues()
	{
		if (!is_copy)
		{
			delete[] q_east;
			delete[] h_east;
			delete[] eta_east;

			delete[] q_west;
			delete[] h_west;
			delete[] eta_west;
		}
	}

} FaceValues;