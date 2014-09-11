#include "gcmethod_2d.h"

int main()
{
	mesh_2d m = mesh_2d("gcmethod_2d.ini");
	m.create_mesh();
	gcmethod_2d test = gcmethod_2d("gcmethod_2d.ini", m);
	test.calculate();
	test.analyze();
	return 0;
}
