CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

PROJECT(fv1-cpp-project)

ADD_EXECUTABLE(fv1-cpp)

TARGET_SOURCES(
    fv1-cpp
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MAIN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input/set_boundary_conditions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input/set_simulation_parameters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input/set_solver_parameters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/add_ghost_cells.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/bed_data_c_property.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/fluxHLL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/friction_implicit.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/friction_update.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/fv1_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/get_bar_values.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/get_dt_CFL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/get_face_values.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/get_modal_values.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/get_nodal_values.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/get_positivity_preserving_nodes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/get_wet_dry_cells.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/h_init_c_property.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator/h_init_overtop.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/output/write_solution_to_file.cpp
)

TARGET_INCLUDE_DIRECTORIES(
    fv1-cpp
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/classes
	${CMAKE_CURRENT_SOURCE_DIR}/input
	${CMAKE_CURRENT_SOURCE_DIR}/operator
	${CMAKE_CURRENT_SOURCE_DIR}/output
)
