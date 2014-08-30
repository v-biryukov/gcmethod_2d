#include "simple2Dgc_method.h"

int main()
{
	simple2Dgcmethod test = simple2Dgcmethod("simple2Dgc_method.ini");
	test.create_mesh();
	test.calculate();
	return 0;
}
