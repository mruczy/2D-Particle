#ifndef V2_H
#define V2_H

#include <math.h>
#include <stdlib.h>

class v2
{
public:
	float x;
	float y;

	v2();
	v2(float xIn, float yIn);
	void randomize();
	virtual ~v2();
};

#endif