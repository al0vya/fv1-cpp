#include "set_error_threshold_epsilon.h"

real set_error_threshold_epsilon()
{
	std::cout << "Please enter an error threshold.\n";

	real user_input_epsilon;

	std::cin >> user_input_epsilon;

	if (!std::cin)
	{
		std::cout << "Error: please rerun and enter a float or double value. Exiting program.\n";

		exit(0);
	}

	return user_input_epsilon;
}