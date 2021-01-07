#include "free_bar_values.h"

void free_bar_values(BarValues& bar_vals)
{
	delete[] bar_vals.h;
	delete[] bar_vals.z;
}