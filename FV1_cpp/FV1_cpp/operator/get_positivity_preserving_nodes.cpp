#include "get_positivity_preserving_nodes.h"

void get_positivity_preserving_nodes
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	FaceValues&           face_vals, 
	StarValues&           star_vals, 
	real*&                delta_west, 
	real*&                delta_east
)
{
	for (int i = 0; i < sim_params.cells + 1; i++)
	{
		real u_west = (face_vals.h_west[i] <= solver_params.tol_dry) ? 0 : face_vals.q_west[i] / face_vals.h_west[i];
		real u_east = (face_vals.h_east[i] <= solver_params.tol_dry) ? 0 : face_vals.q_east[i] / face_vals.h_east[i];

		real z_west = face_vals.eta_west[i] - face_vals.h_west[i];
		real z_east = face_vals.eta_east[i] - face_vals.h_east[i];

		real z_star_intermediate = std::max(z_west, z_east);

		delta_west[i] = std::max(C(0.0), -(face_vals.eta_west[i] - z_star_intermediate));
		delta_east[i] = std::max(C(0.0), -(face_vals.eta_east[i] - z_star_intermediate));

		star_vals.h_west[i] = std::max(C(0.0), face_vals.eta_west[i] - z_star_intermediate);
		star_vals.q_west[i] = u_west * star_vals.h_west[i];

		star_vals.h_east[i] = std::max(C(0.0), face_vals.eta_east[i] - z_star_intermediate);
		star_vals.q_east[i] = u_east * star_vals.h_east[i];

		star_vals.z_west[i] = z_star_intermediate - delta_west[i];
		star_vals.z_east[i] = z_star_intermediate - delta_east[i];
	}
}