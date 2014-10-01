#include "tensor2d.h"
#include <math.h>
#include "mesh_2d/mesh_2d.h"
#include "linear_elastisity_2d.h"



int main()
{
/*
    mesh_2d m = mesh_2d("gcmethod_2d.ini");
    m.create_mesh();

    gcmethod_2d g = gcmethod_2d("gcmethod_2d.ini", &m);
    g.calculate();
    std::ofstream f ("results.txt");
    f << std::fixed << std::setprecision(12) << g.L(1) << '\n' << g.L(2) << '\n' << g.L_inf();
*/
    mesh_2d m = mesh_2d("linear_elastisity.ini");
    m.create_mesh();

    linear_elastisity_2d l = linear_elastisity_2d("linear_elastisity.ini", &m);
    l.calculate();
}

