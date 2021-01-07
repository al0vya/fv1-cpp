#include "set_test_case.h"

int set_test_case()
{
	std::cout << "Please enter a number between 1 and 6 to select a test case.\n"
		         "1: Wet dam break\n"
		         "2: Dry dam break\n"
		         "3: Dry dam break with friction\n"
		         "4: Wet lake-at-rest (C-property)\n"
		         "5: Wet/dry lake-at-rest\n"
		         "6: Building overtopping\n";

	int test_case_selection;

	std::cin >> test_case_selection;

	if (!std::cin || test_case_selection > 6 || test_case_selection < 1)
	{
		std::cout << "Error: please rerun and enter a number between 1 and 6. Exiting program.\n";

		exit(0);
	}

	return test_case_selection;
}