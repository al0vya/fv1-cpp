#include "malloc_face_values.h"

void malloc_face_values
(
	FaceValues&           face_values,
	SimulationParameters& sim_params
)
{
	face_values.q_east   = new real[sim_params.cells + 1];
	face_values.h_east   = new real[sim_params.cells + 1];
	face_values.eta_east = new real[sim_params.cells + 1];

	face_values.q_west   = new real[sim_params.cells + 1];
	face_values.h_west   = new real[sim_params.cells + 1];
	face_values.eta_west = new real[sim_params.cells + 1];
}