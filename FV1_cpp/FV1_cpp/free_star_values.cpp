#include "free_star_values.h"

void free_star_values(StarValues& star_vals)
{
	delete[] star_vals.q_west;
	delete[] star_vals.h_west;

	delete[] star_vals.q_east;
	delete[] star_vals.h_east;

	delete[] star_vals.z_west;
	delete[] star_vals.z_east;
}