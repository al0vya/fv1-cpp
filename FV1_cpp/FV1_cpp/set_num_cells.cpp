#include "set_num_cells.h"

int set_num_cells()
{
	std::cout << "Please select the number of cells.\n";

	int user_input_num_cells;

	std::cin >> user_input_num_cells;

	if (!std::cin || user_input_num_cells < 1)
	{
		std::cout << "Error: please rerun and enter a integer value. Exiting program.\n";

		exit(0);
	}

	return user_input_num_cells;
}