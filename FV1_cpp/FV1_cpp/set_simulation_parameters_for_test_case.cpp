#include "set_simulation_parameters_for_test_case.h"
#include "real.h"

SimulationParameters set_simulation_parameters_for_test_case(int test_case_selection, int number_of_cells)
{
	SimulationParameters simulationParameters;

	switch (test_case_selection)
	{
	case 1:
		simulationParameters.cells = number_of_cells;
		simulationParameters.xmin = 0;
		simulationParameters.xmax = 50;
		simulationParameters.simulationTime = C(2.5);
		simulationParameters.manning = C(0.0);
		break;
	case 2:
		simulationParameters.cells = number_of_cells;
		simulationParameters.xmin = 0;
		simulationParameters.xmax = 50;
		simulationParameters.simulationTime = C(1.3);
		simulationParameters.manning = C(0.0);
		break;
	case 3:
		simulationParameters.cells = number_of_cells;
		simulationParameters.xmin = 0;
		simulationParameters.xmax = 50;
		simulationParameters.simulationTime = C(1.3);
		simulationParameters.manning = C(0.02);
		break;
	case 4:
		simulationParameters.cells = number_of_cells;
		simulationParameters.xmin = 0;
		simulationParameters.xmax = 50;
		simulationParameters.simulationTime = C(0.5);
		simulationParameters.manning = C(0.0);
		break;
	case 5:
		simulationParameters.cells = number_of_cells;
		simulationParameters.xmin = 0;
		simulationParameters.xmax = 50;
		simulationParameters.simulationTime = C(0.5);
		simulationParameters.manning = C(0.0);
		break;
	case 6:
		simulationParameters.cells = number_of_cells;
		simulationParameters.xmin = 0;
		simulationParameters.xmax = 50;
		simulationParameters.simulationTime = C(10.0);
		simulationParameters.manning = C(0.02);
		break;
	default:

		break;
	}

	return simulationParameters;
}