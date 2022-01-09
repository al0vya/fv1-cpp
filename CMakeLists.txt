CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

PROJECT(fv1-cpp-project)

ADD_EXECUTABLE(fv1-cpp)

TARGET_SOURCES(
    fv1-cpp
    PRIVATE
    MAIN.cpp
    input/set_boundary_conditions.cpp
    input/set_simulation_parameters.cpp
    input/set_solver_parameters.cpp
    operator/add_ghost_cells.cpp
    operator/bed_data_c_property.cpp
    operator/fluxHLL.cpp
    operator/friction_implicit.cpp
    operator/friction_update.cpp
    operator/fv1_operator.cpp
    operator/get_bar_values.cpp
    operator/get_dt_CFL.cpp
    operator/get_face_values.cpp
    operator/get_modal_values.cpp
    operator/get_nodal_values.cpp
    operator/get_positivity_preserving_nodes.cpp
    operator/get_wet_dry_cells.cpp
    operator/h_init_c_property.cpp
    operator/h_init_overtop.cpp
    output/write_solution_to_file.cpp
)

TARGET_INCLUDE_DIRECTORIES(
    fv1-cpp
    PRIVATE
    classes
	input
	operator
	output
)
