#include "bed_data_c_property.h"

real bed_data_c_property(real x_int)
{
	real a = x_int;
	real b;

	if (a >= 22 && a < 25)
	{
		b = C(0.05) * a - C(1.1);
	}
	else if (a >= 25 && a <= 28)
	{
		b = C(-0.05) * a + C(1.4);
	}
	else if (a > 8 && a < 12)
	{
		b = C(0.2) - C(0.05) * pow(a - 10, 2);
	}
	else if (a > 39 && a < 46.5)
	{
		b = C(0.3); // whereas here, 0.3 is a double literal, dangerous to cast it to a float
	}
	else
	{
		b = 0; // this is safe because you're casting an int literal to a real
	}

	return b * 10;
}