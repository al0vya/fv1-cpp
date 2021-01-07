#include "free_nodal_values.h"

void free_nodal_values(NodalValues& nodal_vals)
{
	delete[] nodal_vals.q;
	delete[] nodal_vals.h;
	delete[] nodal_vals.z;
}