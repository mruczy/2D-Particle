#include "v2.h"

v2::v2()
{
	randomize();
}

v2::v2(float xIn, float yIn)
{
	x = xIn;
	y = yIn;
}

void v2::randomize()
{
	x = (float)rand() / (float)RAND_MAX;
	y = (float)rand() / (float)RAND_MAX;
}

v2::~v2()
{
}
