///////////////////////////////////////////////////////////
//////
//////	file: convection_equation_solver.cpp
//////  class for solving convection equation in 2d using
//////  grid-characteristic method.
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////

#include "convection_equation_solver.h"

gcmethod_2d::gcmethod_2d()
{
}

gcmethod_2d::gcmethod_2d(std::string path, mesh_2d * mesh_t) : mesh(mesh_t)
{
	read_from_file(path);
	init();
}

gcmethod_2d::gcmethod_2d(struct convection_solver_gc_settings c, mesh_2d * mesh_t) : mesh(mesh_t)
{
    lambda = c.lambda;
    direction = c.direction;
    tau = c.tau;
    number_of_steps = c.number_of_steps;
    N = c.N;
	is_monotonic = c.is_monotonic;
    use_precalculations = c.use_precalculations;
    save_only_main_points = c.save_only_main_points;
    saving_frequency = c.saving_frequency;
	init();
}

void gcmethod_2d::set_number_of_steps(int n)
{
    number_of_steps = n;
}

void gcmethod_2d::save_to_vtk(std::string name)
{
	std::ofstream vtk_file(name.c_str());
	vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
	vtk_file << "DATASET POLYDATA\nPOINTS " << mesh->get_number_of_points() << " float\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
		vtk_file << mesh->get_point(i).x << " " << mesh->get_point(i).y << " "  << 0.0 << "\n";
    if ( !save_only_main_points )
    {
        vtk_file << "\nPOLYGONS " << mesh->get_number_of_elements()*N*N << " " << mesh->get_number_of_elements()*N*N*4 << "\n";
        for (int i = 0; i < mesh->get_number_of_elements(); i++)
        {
            if (N == 1)
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,0) << " " << mesh->get_triangle_point_num(i,1) << " " << mesh->get_triangle_point_num(i,2) << "\n";
            if (N == 2) {
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,0) << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,5) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,1) << " " << mesh->get_triangle_point_num(i,4) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,4) << " " << mesh->get_triangle_point_num(i,2) << " " << mesh->get_triangle_point_num(i,5) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,4) << " " << mesh->get_triangle_point_num(i,5) << "\n";
            }
            if (N == 3) {
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,0) << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,8) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,4) << " " << mesh->get_triangle_point_num(i,9) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,4) << " " << mesh->get_triangle_point_num(i,1) << " " << mesh->get_triangle_point_num(i,5) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,5) << " " << mesh->get_triangle_point_num(i,6) << " " << mesh->get_triangle_point_num(i,9) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,6) << " " << mesh->get_triangle_point_num(i,2) << " " << mesh->get_triangle_point_num(i,7) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,7) << " " << mesh->get_triangle_point_num(i,8) << " " << mesh->get_triangle_point_num(i,9) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,9) << " " << mesh->get_triangle_point_num(i,8) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,4) << " " << mesh->get_triangle_point_num(i,5) << " " << mesh->get_triangle_point_num(i,9) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,6) << " " << mesh->get_triangle_point_num(i,7) << " " << mesh->get_triangle_point_num(i,9) << "\n";
            }
            if (N == 4) {
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,0) << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,11) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,4) << " " << mesh->get_triangle_point_num(i,12) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,4) << " " << mesh->get_triangle_point_num(i,5) << " " << mesh->get_triangle_point_num(i,13) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,5) << " " << mesh->get_triangle_point_num(i,1) << " " << mesh->get_triangle_point_num(i,6) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,6) << " " << mesh->get_triangle_point_num(i,7) << " " << mesh->get_triangle_point_num(i,13) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,7) << " " << mesh->get_triangle_point_num(i,8) << " " << mesh->get_triangle_point_num(i,14) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,8) << " " << mesh->get_triangle_point_num(i,2) << " " << mesh->get_triangle_point_num(i,9) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,9) << " " << mesh->get_triangle_point_num(i,10) << " " << mesh->get_triangle_point_num(i,14) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,10) << " " << mesh->get_triangle_point_num(i,11) << " " << mesh->get_triangle_point_num(i,12) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,12) << " " << mesh->get_triangle_point_num(i,11) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,3) << " " << mesh->get_triangle_point_num(i,13) << " " << mesh->get_triangle_point_num(i,12) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,5) << " " << mesh->get_triangle_point_num(i,6) << " " << mesh->get_triangle_point_num(i,13) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,7) << " " << mesh->get_triangle_point_num(i,14) << " " << mesh->get_triangle_point_num(i,13) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,8) << " " << mesh->get_triangle_point_num(i,9) << " " << mesh->get_triangle_point_num(i,14) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,10) << " " << mesh->get_triangle_point_num(i,12) << " " << mesh->get_triangle_point_num(i,14) << "\n";
                vtk_file << 3 << " " << mesh->get_triangle_point_num(i,12) << " " << mesh->get_triangle_point_num(i,13) << " " << mesh->get_triangle_point_num(i,14) << "\n";
            }
        }
    }
    else
    {
        vtk_file << "\nPOLYGONS " << mesh->get_number_of_elements() << " " << mesh->get_number_of_elements()*4 << "\n";
        for (int i = 0; i < mesh->get_number_of_elements(); i++)
            vtk_file << 3 << " " << mesh->get_triangle_point_num(i,0) << " " << mesh->get_triangle_point_num(i,1) << " " << mesh->get_triangle_point_num(i,2) << "\n";
    }
	vtk_file << "\nPOINT_DATA " << mesh->get_number_of_points() << "\n";;
	vtk_file << "SCALARS z FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
        vtk_file <<  values0.at(i) << "\n";
	}
	vtk_file << "SCALARS exact FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		vtk_file <<  exact_solution(mesh->get_point(i), tau*cur_step) << "\n";
    }
    vtk_file << "SCALARS z_minus_exact FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		vtk_file <<  values0.at(i) - exact_solution(mesh->get_point(i), tau*cur_step) << "\n";
    }
}




void gcmethod_2d::step_any(vector2d step)
{
	vector2d p;
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		p = mesh->get_point(i) + step;
        if (!mesh->is_inside(p)) mesh->make_inside_continuous(p);
        for ( int n : mesh->triangles.at(i) )
            if (mesh->is_inside(p, n))
            {
                values1.at(i) = approximate(p, n);
                break;
            }
	}
}

void gcmethod_2d::fast_step()
{
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        values1[i] = 0;
        for ( int j = 0; j < (N+1)*(N+2)/2; j++ )
            values1[i] += weights[i][j] * values0[mesh->get_triangle_point_num(elements_of_points[i], j)];
    }
    values0.swap(values1);
}

void gcmethod_2d::calculate_weights(vector2d step)
{
	vector2d p;
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		p = mesh->get_point(i) + step;
        if (!mesh->is_inside(p)) mesh->make_inside_continuous(p);
        for ( int n : mesh->triangles.at(i) )
            if (mesh->is_inside(p, n))
            {
                elements_of_points.at(i) = n;
                calculate_weights_of_point(i, p, n);
                break;
            }
	}
}
void gcmethod_2d::calculate_weights_of_point(int pn, vector2d r, int tn)
{
    vector2d ra = mesh->get_triangle_point(tn, 0);
	vector2d rb = mesh->get_triangle_point(tn, 1);
	vector2d rc = mesh->get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	if (N==1)
	{
        weights.at(pn).at(0) = sa;
        weights.at(pn).at(1) = sb;
        weights.at(pn).at(2) = sc;
	}
	else if (N==2)
	{
        weights.at(pn).at(0) = sa*(2*sa-1);
        weights.at(pn).at(1) = sb*(2*sb-1);
        weights.at(pn).at(2) = sc*(2*sc-1);
        weights.at(pn).at(3) = 4*sa*sb;
        weights.at(pn).at(4) = 4*sb*sc;
        weights.at(pn).at(5) = 4*sc*sa;
	}
	else if (N==3)
	{
        weights.at(pn).at(0) = sa*(3*sa-1)*(3*sa-2)/2;
        weights.at(pn).at(1) = sb*(3*sb-1)*(3*sb-2)/2;
        weights.at(pn).at(2) = sc*(3*sc-1)*(3*sc-2)/2;
        weights.at(pn).at(3) = 9*sa*(3*sa-1)*sb/2;
        weights.at(pn).at(4) = 9*sb*(3*sb-1)*sa/2;
        weights.at(pn).at(5) = 9*sb*(3*sb-1)*sc/2;
        weights.at(pn).at(6) = 9*sc*(3*sc-1)*sb/2;
        weights.at(pn).at(7) = 9*sc*(3*sc-1)*sa/2;
        weights.at(pn).at(8) = 9*sa*(3*sa-1)*sc/2;
        weights.at(pn).at(9) = 27*sa*sb*sc;
	}
	else if (N==4)
	{
        weights.at(pn).at(0) = sa*(4*sa-1)*(2*sa-1)*(4*sa-3)/3;
        weights.at(pn).at(1) = sb*(4*sb-1)*(2*sb-1)*(4*sb-3)/3;
        weights.at(pn).at(2) = sc*(4*sc-1)*(2*sc-1)*(4*sc-3)/3;
        weights.at(pn).at(3) = 16*sa*(4*sa-1)*(2*sa-1)*sb/3;
        weights.at(pn).at(4) = 4*sa*(4*sa-1)*sb*(4*sb-1);
        weights.at(pn).at(5) = 16*sb*(4*sb-1)*(2*sb-1)*sa/3;
        weights.at(pn).at(6) = 16*sb*(4*sb-1)*(2*sb-1)*sc/3;
        weights.at(pn).at(7) = 4*sb*(4*sb-1)*sc*(4*sc-1);
        weights.at(pn).at(8) = 16*sc*(4*sc-1)*(2*sc-1)*sb/3;
        weights.at(pn).at(9) = 16*sc*(4*sc-1)*(2*sc-1)*sa/3;
        weights.at(pn).at(10) = 4*sa*(4*sa-1)*sc*(4*sc-1);
        weights.at(pn).at(11) = 16*sa*(4*sa-1)*(2*sa-1)*sc/3;
        weights.at(pn).at(12) = 32*sa*(4*sa-1)*sb*sc;
        weights.at(pn).at(13) = 32*sb*(4*sb-1)*sc*sa;
        weights.at(pn).at(14) = 32*sc*(4*sc-1)*sa*sb;
	}
}

bool gcmethod_2d::min_max_check(double z, int tn)
{
    return   (  z > values0.at(mesh->get_triangle_point_num(tn,0)) + eps_z
             && z > values0.at(mesh->get_triangle_point_num(tn,1)) + eps_z
             && z > values0.at(mesh->get_triangle_point_num(tn,2)) + eps_z )
           ||(  z < values0.at(mesh->get_triangle_point_num(tn,0)) - eps_z
             && z < values0.at(mesh->get_triangle_point_num(tn,1)) - eps_z
             && z < values0.at(mesh->get_triangle_point_num(tn,2)) - eps_z );
}

double gcmethod_2d::approximate(vector2d p, int tn)
{
    double result;
	switch (N)
	{
		case 1: result = approximate_linear(p, tn); break;
		case 2: result = approximate_quadratic_special(p, tn); break;
		case 3: result = approximate_cubic(p, tn); break;
		case 4: result = approximate_quartic(p, tn); break;
		default: break;
	}
	if (is_monotonic && min_max_check(result, tn))
        return approximate_linear(p, tn);
    else
        return result;
}

double gcmethod_2d::approximate_linear(vector2d r, int tn)
{
	vector2d ra = mesh->get_triangle_point(tn, 0);
	vector2d rb = mesh->get_triangle_point(tn, 1);
	vector2d rc = mesh->get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	double v = 0;
	v += sa * values0.at(mesh->get_triangle_point_num(tn, 0));
	v += sb * values0.at(mesh->get_triangle_point_num(tn, 1));
	v += sc * values0.at(mesh->get_triangle_point_num(tn, 2));
	return v;
}

double gcmethod_2d::approximate_quadratic(vector2d r, int tn)
{
    vector2d ra = mesh->get_triangle_point(tn, 0);
	vector2d rb = mesh->get_triangle_point(tn, 1);
	vector2d rc = mesh->get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	double v = 0;
	v += sa*(2*sa-1) * values0.at(mesh->get_triangle_point_num(tn, 0));
	v += sb*(2*sb-1) * values0.at(mesh->get_triangle_point_num(tn, 1));
	v += sc*(2*sc-1) * values0.at(mesh->get_triangle_point_num(tn, 2));
	v += 4*sa*sb * values0.at(mesh->get_triangle_point_num(tn, 3));
	v += 4*sb*sc * values0.at(mesh->get_triangle_point_num(tn, 4));
	v += 4*sc*sa * values0.at(mesh->get_triangle_point_num(tn, 5));//std::cout << v << "\n";
	return v;
}

double gcmethod_2d::approximate_cubic(vector2d r, int tn)
{
    vector2d ra = mesh->get_triangle_point(tn, 0);
	vector2d rb = mesh->get_triangle_point(tn, 1);
	vector2d rc = mesh->get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	double v = 0;
	v += sa*(3*sa-1)*(3*sa-2)/2 * values0.at(mesh->get_triangle_point_num(tn, 0));
	v += sb*(3*sb-1)*(3*sb-2)/2 * values0.at(mesh->get_triangle_point_num(tn, 1));
	v += sc*(3*sc-1)*(3*sc-2)/2 * values0.at(mesh->get_triangle_point_num(tn, 2));
	v += 9*sa*(3*sa-1)*sb/2     * values0.at(mesh->get_triangle_point_num(tn, 3));
	v += 9*sb*(3*sb-1)*sa/2     * values0.at(mesh->get_triangle_point_num(tn, 4));
    v += 9*sb*(3*sb-1)*sc/2     * values0.at(mesh->get_triangle_point_num(tn, 5));
    v += 9*sc*(3*sc-1)*sb/2     * values0.at(mesh->get_triangle_point_num(tn, 6));
	v += 9*sc*(3*sc-1)*sa/2     * values0.at(mesh->get_triangle_point_num(tn, 7));
	v += 9*sa*(3*sa-1)*sc/2     * values0.at(mesh->get_triangle_point_num(tn, 8));
	v += 27*sa*sb*sc            * values0.at(mesh->get_triangle_point_num(tn, 9));
	return v;
}

double gcmethod_2d::approximate_quartic(vector2d r, int tn)
{
    vector2d ra = mesh->get_triangle_point(tn, 0);
	vector2d rb = mesh->get_triangle_point(tn, 1);
	vector2d rc = mesh->get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;

	double v = 0;
	v += sa*(4*sa-1)*(2*sa-1)*(4*sa-3)/3 * values0.at(mesh->get_triangle_point_num(tn, 0));
	v += sb*(4*sb-1)*(2*sb-1)*(4*sb-3)/3 * values0.at(mesh->get_triangle_point_num(tn, 1));
	v += sc*(4*sc-1)*(2*sc-1)*(4*sc-3)/3 * values0.at(mesh->get_triangle_point_num(tn, 2));
	v += 16*sa*(4*sa-1)*(2*sa-1)*sb/3 * values0.at(mesh->get_triangle_point_num(tn, 3));
	v += 4*sa*(4*sa-1)*sb*(4*sb-1)    * values0.at(mesh->get_triangle_point_num(tn, 4));
	v += 16*sb*(4*sb-1)*(2*sb-1)*sa/3 * values0.at(mesh->get_triangle_point_num(tn, 5));
	v += 16*sb*(4*sb-1)*(2*sb-1)*sc/3 * values0.at(mesh->get_triangle_point_num(tn, 6));
	v += 4*sb*(4*sb-1)*sc*(4*sc-1)    * values0.at(mesh->get_triangle_point_num(tn, 7));
    v += 16*sc*(4*sc-1)*(2*sc-1)*sb/3 * values0.at(mesh->get_triangle_point_num(tn, 8));
    v += 16*sc*(4*sc-1)*(2*sc-1)*sa/3 * values0.at(mesh->get_triangle_point_num(tn, 9));
    v += 4*sa*(4*sa-1)*sc*(4*sc-1)    * values0.at(mesh->get_triangle_point_num(tn, 10));
	v += 16*sa*(4*sa-1)*(2*sa-1)*sc/3 * values0.at(mesh->get_triangle_point_num(tn, 11));
	v += 32*sa*(4*sa-1)*sb*sc         * values0.at(mesh->get_triangle_point_num(tn, 12));
    v += 32*sb*(4*sb-1)*sc*sa         * values0.at(mesh->get_triangle_point_num(tn, 13));
    v += 32*sc*(4*sc-1)*sa*sb         * values0.at(mesh->get_triangle_point_num(tn, 14));
	return v;
}

void gcmethod_2d::aproximate_quadratic_special_small(std::vector<double> & zs, int tn, vector2d p1, vector2d p2, vector2d p3, double z1, double z2, double z3)
{
    zs.push_back(z1);
    zs.push_back(z2);
    zs.push_back(z3);
    zs.push_back(approximate_quadratic((p1+p2)/2, tn));
    zs.push_back(approximate_quadratic((p1+p3)/2, tn));
    zs.push_back(approximate_quadratic((p3+p2)/2, tn));

    if (  (zs.at(3) > z1 && zs.at(3) > z2)
        ||(zs.at(3) < z1 && zs.at(3) < z2)  )
        zs.at(3) = (z1+z2)/2;
    if (  (zs.at(4) > z1 && zs.at(4) > z3)
        ||(zs.at(4) < z1 && zs.at(4) < z3)  )
        zs.at(4) = (z3+z1)/2;
    if (  (zs.at(5) > z3 && zs.at(5) > z2)
        ||(zs.at(5) < z3 && zs.at(5) < z2)  )
        zs.at(5) = (z3+z2)/2;
}

double gcmethod_2d::approximate_quadratic_special(vector2d r, int tn)
{
    vector2d ra = mesh->get_triangle_point(tn, 0);
	vector2d rb = mesh->get_triangle_point(tn, 1);
	vector2d rc = mesh->get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	vector2d p1, p2, p3;
	std::vector<double> z;
    if (sa >= 0.5)
    {
        p1 = ra;
        p2 = ra + (rb-ra)/2;
        p3 = ra + (rc-ra)/2;
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, values0.at(mesh->elements.at(tn).at(0)),
                                                              values0.at(mesh->elements.at(tn).at(3)),
                                                              values0.at(mesh->elements.at(tn).at(5)));
    }
    else if (sb >= 0.5)
    {
        p1 = rb;
        p2 = rb + (rc-rb)/2;
        p3 = rb + (ra-rb)/2;
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, values0.at(mesh->elements.at(tn).at(1)),
                                                              values0.at(mesh->elements.at(tn).at(4)),
                                                              values0.at(mesh->elements.at(tn).at(3)));
    }
    else if (sc >= 0.5)
    {
        p1 = rc;
        p2 = rc + (ra-rc)/2;
        p3 = rc + (rb-rc)/2;
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, values0.at(mesh->elements.at(tn).at(2)),
                                                              values0.at(mesh->elements.at(tn).at(5)),
                                                              values0.at(mesh->elements.at(tn).at(4)));
    }
    else
    {
        p1 = (ra+rb)/2;
        p2 = (rb+rc)/2;
        p3 = (rc+ra)/2;
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, values0.at(mesh->elements.at(tn).at(3)),
                                                              values0.at(mesh->elements.at(tn).at(4)),
                                                              values0.at(mesh->elements.at(tn).at(5)));
    }
	double spa = Vec(p3 - p2, r - p2)/2;
	double spb = Vec(p1 - p3, r - p3)/2;
	double spc = Vec(p2 - p1, r - p1)/2;
	double sp = spa + spb + spc;
	spa /= sp; spb /= sp; spc /= sp;
	double v = 0;
	v += spa*(2*spa-1) * z.at(0);
	v += spb*(2*spb-1) * z.at(1);
	v += spc*(2*spc-1) * z.at(2);
	v += 4*spa*spb * z.at(3);
	v += 4*spc*spa * z.at(4);
	v += 4*spb*spc * z.at(5);
	return v;
}

void gcmethod_2d::step()
{
	step_any(-tau * lambda * direction);
	values1.swap(values0);
	cur_step++;
}



void gcmethod_2d::calculate()
{
	if ( use_precalculations )
        calculate_weights(-tau * lambda * direction);
	for (int i = 0; i < number_of_steps; i++)
	{
        //DrawProgressBar(40, i*1.0/number_of_steps);
        if (saving_frequency > 0 && i % saving_frequency == 0 )
            save_to_vtk("out/out_" + std::to_string(i) + ".vtk");
		if ( use_precalculations )
            fast_step();
        else
            step();
	}
	//DrawProgressBar(40, 1.0);
	if (saving_frequency > 0)
        save_to_vtk("out/out_" + std::to_string(number_of_steps) + ".vtk");
}

void gcmethod_2d::init()
{
    cur_step = 0;
	values0.resize(mesh->get_number_of_points());
	values1.resize(mesh->get_number_of_points());
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
		values0.at(i) = initial_conditions(mesh->get_point(i));
    if ( use_precalculations )
    {
        weights.resize(mesh->get_number_of_points(), std::vector<double>((N+1)*(N+2)/2));
        elements_of_points.resize(mesh->get_number_of_points());
    }
}



double gcmethod_2d::exact_solution( vector2d p, double t )
{
    vector2d r = p - t * lambda * direction;
    while (!mesh->is_inside(r))
        mesh->make_inside_continuous(r);
    return initial_conditions(r);
}

double gcmethod_2d::L_inf()
{
    double max_err = fabs(values0.at(0) - exact_solution(mesh->get_point(0), tau*number_of_steps));
    for ( int i = 1; i < mesh->get_number_of_points(); i++ )
    {
        if ( max_err < fabs(values0.at(i) - exact_solution(mesh->get_point(i), tau*number_of_steps)) )
            max_err = fabs(values0.at(i) - exact_solution(mesh->get_point(i), tau*number_of_steps));
    }
    return max_err;
}

double gcmethod_2d::L(int n)
{
    double sum = 0, temp;
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        temp = pow(fabs(values0.at(i) - exact_solution(mesh->get_point(i), tau*number_of_steps)), n);
        if ( mesh->get_is_structured() )
        {
            temp *= mesh->get_h() * mesh->get_h();
            if (mesh->get_point(i).x < -mesh->get_size_x()/2 + eps_xy || mesh->get_point(i).x > mesh->get_size_x()/2 - eps_xy)
                temp /= 2;
            if (mesh->get_point(i).y < -mesh->get_size_y()/2 + eps_xy || mesh->get_point(i).y > mesh->get_size_y()/2 - eps_xy)
                temp /= 2;
            sum +=  temp;
        }
        else
        {
            sum += mesh->voronoi_areas.at(i) * temp;
        }
    }
    return pow(sum, 1.0/n) / mesh->get_size_x() / mesh->get_size_y();
}

void gcmethod_2d::analyze()
{
    std::cout << "Results:" << std::endl;
    std::cout << "Order of the approximation is " << N << std::endl;
    std::cout << "L_infinity is equal to " << std::fixed << std::setprecision(10) << L_inf() << std::endl;
    std::cout << "L_1 is equal to " << std::fixed << std::setprecision(10) << L(1) << std::endl;
    std::cout << "L_2 is equal to " << std::fixed << std::setprecision(10) << L(2) << std::endl;
}

void gcmethod_2d::read_from_file(std::string path)
{
    using std::string;
    using std::cout;
    using std::cin;
    using std::endl;
    boost::property_tree::ptree pt;
    try
    {
        boost::property_tree::read_ini(path, pt);
    }
    catch (boost::property_tree::ini_parser_error& error)
    {
        cout
            << error.message() << ": "
            << error.filename() << ", line "
            << error.line() << endl;
        cout << "Error! Press any key to close." << endl;
        std::cin.get();
        std::exit(1);
    }
    lambda = pt.get<REAL>("Equation.lambda");
    std::string direction_str = pt.get<std::string>("Equation.direction");
    std::stringstream ss;
    ss << direction_str;
    ss >> direction.x >> direction.y;
    tau = pt.get<REAL>("Method.tau");
    number_of_steps = pt.get<int>("Method.number_of_steps");
    N = pt.get<int>("Method.N");
    std::string is_monotonic_str = pt.get<std::string>("Method.is_monotonic");
    is_monotonic = (is_monotonic_str == "true" || is_monotonic_str == "TRUE" || is_monotonic_str == "True");
    std::string use_precalculations_str = pt.get<std::string>("Method.use_precalculations");
    use_precalculations = (use_precalculations_str == "true" || use_precalculations_str == "TRUE" || use_precalculations_str == "True");
    std::string save_only_main_points_str = pt.get<std::string>("Method.save_only_main_points");
    save_only_main_points = (save_only_main_points_str == "true" || save_only_main_points_str == "TRUE" || save_only_main_points_str == "True");
    saving_frequency = pt.get<int>("Method.saving_frequency");
}

