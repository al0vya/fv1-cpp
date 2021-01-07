#include "free_face_values.h"

void free_face_values(FaceValues& face_vals)
{
	delete[] face_vals.q_east;
	delete[] face_vals.h_east;
	delete[] face_vals.eta_east;

	delete[] face_vals.q_west;
	delete[] face_vals.h_west;
	delete[] face_vals.eta_west;
}