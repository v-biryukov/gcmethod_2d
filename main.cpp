#include "gcmethod_2d.h"

int main()
{
	gcmethod_2d test = gcmethod_2d("gcmethod_2d.ini");
	test.create_mesh();
	test.calculate();
	return 0;
}
