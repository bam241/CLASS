#ifndef _CLASSMETHOD_
#define _CLASSMETHOD_

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>

struct my_tolower
{
	char operator()(char c) const
	{
		return std::tolower(static_cast<unsigned char>(c));
	}
};

//To Lower Case, convert any string in lower case
static std::string tlc(string data)
{
	transform(data.begin(), data.end(), data.begin(), my_tolower());
	return data;
};



static float random(float a, float b) //peak random numebr between a and b
{
	float range = pow(2., 31);
	srand(time(NULL)); //initialize the srand
	return (float)a + (float)(b-a)*rand()/range;
};

static std::string dtoa(double num)
{
	std::ostringstream os(std::ostringstream::out);
	os << setprecision(3) << num;
	return os.str();
};


#endif
