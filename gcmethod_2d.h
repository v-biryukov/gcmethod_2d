///////////////////////////////////////////////////////////
//////
//////  class for solving convection equation in 2d using
//////  grid-characteristic method.
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////


#pragma once

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#define ANSI_DECLARATORS

extern "C"
{
#include "triangle.h"
}

#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iomanip>
#include "vector2.h"
#include "mesh_2d.h"



class gcmethod_2d
{
	std::string path;
	double lambda_x, lambda_y;
	double tau;
	int number_of_steps;
	int cur_step;

	double eps_xy = 1e-10;
	double eps_z = 1e-10;
	mesh_2d & mesh;
	std::vector<double> values0;
	std::vector<double> values1;

	std::vector<std::vector<double> > additional_z0;
	std::vector<std::vector<double> > additional_z1;

	std::vector<std::vector<double> > weights;

public:
    int N;
	bool is_monotonic;
	gcmethod_2d(std::string path, mesh_2d & mesh_t);
	void calculate();
	void save_to_vtk(std::string name);
	void analyze();
    double L_inf();
    double L(int n);
private:
	void init();
	vector2d get_additional_point(int n, int k);
	void read_from_file(std::string path);
	double initial_conditions(double x, double y)
        {if (x*x + y*y < 9) return 2.5; else return 0.0;};//return 2*exp(-(x*x+y*y)/5);};//{if (x*x + y*y < 5) return 2.5; else return 0.0;};//return 2.0/(1 + x*x + y*y);};
	double initial_conditions(vector2d v) {return initial_conditions(v.x, v.y);};
	void step_any(vector2d step);
	double approximate(vector2d p, int tn);
    double prepare_approximate_linear(vector2d p, int tn);
	double approximate_linear(vector2d p, int tn);
	double approximate_quadratic(vector2d p, int tn);
    double approximate_cubic(vector2d p, int tn);
    double approximate_quartic(vector2d r, int tn);
    double approximate_quadratic_special(vector2d r, int tn);
    void aproximate_quadratic_special_small(std::vector<double> & zs,  int tn, vector2d p1, vector2d p2, vector2d p3, double z1, double z2, double z3);
    double approximate_any(vector2d p, int tn);
    bool min_max_check(double z, int tn);
    double exact_solution(vector2d p, double t);

	void step();
};

gcmethod_2d::gcmethod_2d(std::string path, mesh_2d & mesh_t) : mesh(mesh_t)
{
	read_from_file(path);
	cur_step = 0;

}

void gcmethod_2d::save_to_vtk(std::string name)
{
	std::ofstream vtk_file(name.c_str());
	vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
	vtk_file << "DATASET POLYDATA\nPOINTS " << mesh.get_number_of_points() << " float\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
		vtk_file << mesh.get_point(i).x << " " << mesh.get_point(i).y << " "  << 0.0 << "\n";
	vtk_file << "\nPOLYGONS " << mesh.get_number_of_triangles()*N*N << " " << mesh.get_number_of_triangles()*N*N*4 << "\n";
	for (int i = 0; i < mesh.get_number_of_triangles(); i++)
	{
		//vtk_file << 3 << " " << mesh.get_triangle_point_num(i,0) << " " << mesh.get_triangle_point_num(i,1) << " " << mesh.get_triangle_point_num(i,2) << "\n";
        if (N == 1)
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,0) << " " << mesh.get_triangle_point_num(i,1) << " " << mesh.get_triangle_point_num(i,2) << "\n";
        if (N == 2) {
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,0) << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,5) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,1) << " " << mesh.get_triangle_point_num(i,4) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,4) << " " << mesh.get_triangle_point_num(i,2) << " " << mesh.get_triangle_point_num(i,5) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,4) << " " << mesh.get_triangle_point_num(i,5) << "\n";
        }
        if (N == 3) {
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,0) << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,8) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,4) << " " << mesh.get_triangle_point_num(i,9) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,4) << " " << mesh.get_triangle_point_num(i,1) << " " << mesh.get_triangle_point_num(i,5) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,5) << " " << mesh.get_triangle_point_num(i,6) << " " << mesh.get_triangle_point_num(i,9) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,6) << " " << mesh.get_triangle_point_num(i,2) << " " << mesh.get_triangle_point_num(i,7) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,7) << " " << mesh.get_triangle_point_num(i,8) << " " << mesh.get_triangle_point_num(i,9) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,9) << " " << mesh.get_triangle_point_num(i,8) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,4) << " " << mesh.get_triangle_point_num(i,5) << " " << mesh.get_triangle_point_num(i,9) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,6) << " " << mesh.get_triangle_point_num(i,7) << " " << mesh.get_triangle_point_num(i,9) << "\n";
        }
        if (N == 4) {
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,0) << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,11) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,4) << " " << mesh.get_triangle_point_num(i,12) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,4) << " " << mesh.get_triangle_point_num(i,5) << " " << mesh.get_triangle_point_num(i,13) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,5) << " " << mesh.get_triangle_point_num(i,1) << " " << mesh.get_triangle_point_num(i,6) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,6) << " " << mesh.get_triangle_point_num(i,7) << " " << mesh.get_triangle_point_num(i,13) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,7) << " " << mesh.get_triangle_point_num(i,8) << " " << mesh.get_triangle_point_num(i,14) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,8) << " " << mesh.get_triangle_point_num(i,2) << " " << mesh.get_triangle_point_num(i,9) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,9) << " " << mesh.get_triangle_point_num(i,10) << " " << mesh.get_triangle_point_num(i,14) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,10) << " " << mesh.get_triangle_point_num(i,11) << " " << mesh.get_triangle_point_num(i,12) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,12) << " " << mesh.get_triangle_point_num(i,11) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,3) << " " << mesh.get_triangle_point_num(i,13) << " " << mesh.get_triangle_point_num(i,12) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,5) << " " << mesh.get_triangle_point_num(i,6) << " " << mesh.get_triangle_point_num(i,13) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,7) << " " << mesh.get_triangle_point_num(i,14) << " " << mesh.get_triangle_point_num(i,13) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,8) << " " << mesh.get_triangle_point_num(i,9) << " " << mesh.get_triangle_point_num(i,14) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,10) << " " << mesh.get_triangle_point_num(i,12) << " " << mesh.get_triangle_point_num(i,14) << "\n";
            vtk_file << 3 << " " << mesh.get_triangle_point_num(i,12) << " " << mesh.get_triangle_point_num(i,13) << " " << mesh.get_triangle_point_num(i,14) << "\n";
        }
    }
	vtk_file << "\nPOINT_DATA " << mesh.get_number_of_points() << "\n";;
	vtk_file << "SCALARS z FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{
        vtk_file <<  values0.at(i) << "\n";
	}
	vtk_file << "SCALARS exact FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{
		vtk_file <<  exact_solution(mesh.get_point(i), tau*cur_step) << "\n";
    }
    vtk_file << "SCALARS z_minus_exact FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{
		vtk_file <<  values0.at(i) - exact_solution(mesh.get_point(i), tau*cur_step) << "\n";
    }
}




void gcmethod_2d::step_any(vector2d step)
{
	vector2d p;
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{
		p = mesh.get_point(i) + step;
		if ((p.x-0)*(p.x-0) + p.y*p.y < 1)
            int a = 0;
		if (!mesh.is_inside(p)) mesh.make_inside_vector(p);
        for ( int n : mesh.triangles.at(i) )
            if (mesh.is_inside(p, n))
            {
                values1.at(i) = approximate(p, n);
                break;
            }
	}
}

bool gcmethod_2d::min_max_check(double z, int tn)
{
    return   (  z > values0.at(mesh.get_triangle_point_num(tn,0)) + eps_z
             && z > values0.at(mesh.get_triangle_point_num(tn,1)) + eps_z
             && z > values0.at(mesh.get_triangle_point_num(tn,2)) + eps_z )
           ||(  z < values0.at(mesh.get_triangle_point_num(tn,0)) - eps_z
             && z < values0.at(mesh.get_triangle_point_num(tn,1)) - eps_z
             && z < values0.at(mesh.get_triangle_point_num(tn,2)) - eps_z );
}

double gcmethod_2d::approximate(vector2d p, int tn)
{
    double result;
	switch (N)
	{
		case 1: result = prepare_approximate_linear(p, tn); break;
		case 2: result = approximate_quadratic(p, tn); break;
		case 3: result = approximate_cubic(p, tn); break;
		case 4: result = approximate_quartic(p, tn); break;
		default: break;
	}
	if (is_monotonic && min_max_check(result, tn))
        return prepare_approximate_linear(p, tn);
    else
        return result;
}


double gcmethod_2d::prepare_approximate_linear(vector2d p, int tn)
{
	double x1 = mesh.get_triangle_point(tn, 0).x; double y1 = mesh.get_triangle_point(tn, 0).y; double z1 = values0.at(mesh.get_triangle_point_num(tn, 0));
	double x2 = mesh.get_triangle_point(tn, 1).x; double y2 = mesh.get_triangle_point(tn, 1).y; double z2 = values0.at(mesh.get_triangle_point_num(tn, 1));
	double x3 = mesh.get_triangle_point(tn, 2).x; double y3 = mesh.get_triangle_point(tn, 2).y; double z3 = values0.at(mesh.get_triangle_point_num(tn, 2));

	//double vx = (y3-y1)*(z3-z2) - (z3-z1)*(y3-y2);
	//double vy = (z3-z1)*(x3-x2) - (x3-x1)*(z3-z2);
	double vz = (x3-x1)*(y3-y2) - (y3-y1)*(x3-x2);

    //weights.at(tn).at(0) = ((x1-p.x)*(y3-y2) + (y1-p.y)*(x2-x3) + vz) / vz;
    //weights.at(tn).at(1) = ((x1-p.x)*(y1-y3) + (y1-p.y)*(x3-x1)) / vz;
    //weights.at(tn).at(2) = ((x1-p.x)*(y2-y1) + (y1-p.y)*(x1-x2)) / vz;
    return ((x1-p.x)*(y3-y2) + (y1-p.y)*(x2-x3) + vz) / vz*z1 + ((x1-p.x)*(y1-y3) + (y1-p.y)*(x3-x1)) / vz*z2 + ((x1-p.x)*(y2-y1) + (y1-p.y)*(x1-x2)) / vz*z3;
}

double gcmethod_2d::approximate_linear(vector2d p, int tn)
{
	double z1 = values0.at(mesh.get_triangle_point_num(tn, 0));
	double z2 = values0.at(mesh.get_triangle_point_num(tn, 1));
	double z3 = values0.at(mesh.get_triangle_point_num(tn, 2));
	return weights.at(tn).at(0)*z1 + weights.at(tn).at(1)*z2 + weights.at(tn).at(2)*z3;
}

double gcmethod_2d::approximate_quadratic(vector2d r, int tn)
{
    vector2d ra = mesh.get_triangle_point(tn, 0);
	vector2d rb = mesh.get_triangle_point(tn, 1);
	vector2d rc = mesh.get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	double v = 0;
	v += sa*(2*sa-1) * values0.at(mesh.get_triangle_point_num(tn, 0));
	v += sb*(2*sb-1) * values0.at(mesh.get_triangle_point_num(tn, 1));
	v += sc*(2*sc-1) * values0.at(mesh.get_triangle_point_num(tn, 2));
	v += 4*sa*sb * values0.at(mesh.get_triangle_point_num(tn, 3));
	v += 4*sb*sc * values0.at(mesh.get_triangle_point_num(tn, 4));
	v += 4*sc*sa * values0.at(mesh.get_triangle_point_num(tn, 5));//std::cout << v << "\n";
    if (r.x*r.x + r.y*r.y < 6 && v < 0.5)
        int b = 1;
	return v;
}

double gcmethod_2d::approximate_cubic(vector2d r, int tn)
{
    vector2d ra = mesh.get_triangle_point(tn, 0);
	vector2d rb = mesh.get_triangle_point(tn, 1);
	vector2d rc = mesh.get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	double v = 0;
	v += sa*(3*sa-1)*(3*sa-2)/2 * values0.at(mesh.get_triangle_point_num(tn, 0));
	v += sb*(3*sb-1)*(3*sb-2)/2 * values0.at(mesh.get_triangle_point_num(tn, 1));
	v += sc*(3*sc-1)*(3*sc-2)/2 * values0.at(mesh.get_triangle_point_num(tn, 2));
	v += 9*sa*(3*sa-1)*sb/2     * values0.at(mesh.get_triangle_point_num(tn, 3));
	v += 9*sb*(3*sb-1)*sa/2     * values0.at(mesh.get_triangle_point_num(tn, 4));
    v += 9*sb*(3*sb-1)*sc/2     * values0.at(mesh.get_triangle_point_num(tn, 5));
    v += 9*sc*(3*sc-1)*sb/2     * values0.at(mesh.get_triangle_point_num(tn, 6));
	v += 9*sc*(3*sc-1)*sa/2     * values0.at(mesh.get_triangle_point_num(tn, 7));
	v += 9*sa*(3*sa-1)*sc/2     * values0.at(mesh.get_triangle_point_num(tn, 8));
	v += 27*sa*sb*sc            * values0.at(mesh.get_triangle_point_num(tn, 9));
	return v;
}

double gcmethod_2d::approximate_quartic(vector2d r, int tn)
{
    vector2d ra = mesh.get_triangle_point(tn, 0);
	vector2d rb = mesh.get_triangle_point(tn, 1);
	vector2d rc = mesh.get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;

	double v = 0;
	v += sa*(4*sa-1)*(2*sa-1)*(4*sa-3)/3 * values0.at(mesh.get_triangle_point_num(tn, 0));
	v += sb*(4*sb-1)*(2*sb-1)*(4*sb-3)/3 * values0.at(mesh.get_triangle_point_num(tn, 1));
	v += sc*(4*sc-1)*(2*sc-1)*(4*sc-3)/3 * values0.at(mesh.get_triangle_point_num(tn, 2));
	v += 16*sa*(4*sa-1)*(2*sa-1)*sb/3 * values0.at(mesh.get_triangle_point_num(tn, 3));
	v += 4*sa*(4*sa-1)*sb*(4*sb-1)    * values0.at(mesh.get_triangle_point_num(tn, 4));
	v += 16*sb*(4*sb-1)*(2*sb-1)*sa/3 * values0.at(mesh.get_triangle_point_num(tn, 5));
	v += 16*sb*(4*sb-1)*(2*sb-1)*sc/3 * values0.at(mesh.get_triangle_point_num(tn, 6));
	v += 4*sb*(4*sb-1)*sc*(4*sc-1)    * values0.at(mesh.get_triangle_point_num(tn, 7));
    v += 16*sc*(4*sc-1)*(2*sc-1)*sb/3 * values0.at(mesh.get_triangle_point_num(tn, 8));
    v += 16*sc*(4*sc-1)*(2*sc-1)*sa/3 * values0.at(mesh.get_triangle_point_num(tn, 9));
    v += 4*sa*(4*sa-1)*sc*(4*sc-1)    * values0.at(mesh.get_triangle_point_num(tn, 10));
	v += 16*sa*(4*sa-1)*(2*sa-1)*sc/3 * values0.at(mesh.get_triangle_point_num(tn, 11));
	v += 32*sa*(4*sa-1)*sb*sc         * values0.at(mesh.get_triangle_point_num(tn, 12));
    v += 32*sb*(4*sb-1)*sc*sa         * values0.at(mesh.get_triangle_point_num(tn, 13));
    v += 32*sc*(4*sc-1)*sa*sb         * values0.at(mesh.get_triangle_point_num(tn, 14));
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
    vector2d ra = mesh.get_triangle_point(tn, 0);
	vector2d rb = mesh.get_triangle_point(tn, 1);
	vector2d rc = mesh.get_triangle_point(tn, 2);
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
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, values0.at(mesh.get_triangle_point_num(tn, 0)), additional_z0.at(tn).at(0), additional_z0.at(tn).at(1));
    }
    else if (sb >= 0.5)
    {
        p1 = rb;
        p2 = rb + (rc-rb)/2;
        p3 = rb + (ra-rb)/2;
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, values0.at(mesh.get_triangle_point_num(tn, 1)), additional_z0.at(tn).at(2), additional_z0.at(tn).at(0));
    }
    else if (sc >= 0.5)
    {
        p1 = rc;
        p2 = rc + (ra-rc)/2;
        p3 = rc + (rb-rc)/2;
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, values0.at(mesh.get_triangle_point_num(tn, 2)), additional_z0.at(tn).at(1), additional_z0.at(tn).at(2));
    }
    else
    {
        p1 = (ra+rb)/2;
        p2 = (rb+rc)/2;
        p3 = (rc+ra)/2;
        aproximate_quadratic_special_small(z, tn, p1, p2, p3, additional_z0.at(tn).at(0), additional_z0.at(tn).at(2), additional_z0.at(tn).at(1));
    }
	double spa = Vec(p3 - p2, r - p2)/2;
	double spb = Vec(p1 - p3, r - p3)/2;
	double spc = Vec(p2 - p1, r - p1)/2;
	double sp = spa + spb + spc;
	spa /= sp; spb /= sp; spc /= sp;
    //if (r.x*r.x + r.y*r.y < 1)
    //    std::cout << z.at(0) << " " << z.at(1) << " " << z.at(2) << " " << z.at(3) << " " << z.at(4) << " " << z.at(5) << "\n";
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
	// time step for the equation u_t + lx*u_x + ly*u_y = 0;
	step_any(vector2d(-tau * lambda_x, 0));
	values1.swap(values0);
	additional_z1.swap(additional_z0);
	step_any(vector2d(0, -tau * lambda_y));
	values1.swap(values0);
	additional_z1.swap(additional_z0);
	cur_step++;
}


void DrawProgressBar(int len, double percent) {
  std::cout << "\x1B[2K"; // Erase the entire current line.
  std::cout << "\x1B[0E"; // Move to the beginning of the current line.
  std::string progress;
  for (int i = 0; i < len; ++i) {
    if (i < static_cast<int>(len * percent)) {
      progress += "=";
    } else {
      progress += " ";
    }
  }
  std::cout << "[" << progress << "] " << (static_cast<int>(100 * percent)) << "%";
  flush(std::cout); // Required.
}

void gcmethod_2d::calculate()
{
	init();
	for (int i = 0; i < number_of_steps; i++)
	{
        // saving every step
        DrawProgressBar(40, i*1.0/number_of_steps);
		save_to_vtk("out/out_" + std::to_string(i) + ".vtk");
		step();

	}
	DrawProgressBar(40, 1.0);
	save_to_vtk("out/out_" + std::to_string(number_of_steps) + ".vtk");
}

void gcmethod_2d::init()
{
	values0.resize(mesh.get_number_of_points());
	values1.resize(mesh.get_number_of_points());
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
		values0.at(i) = initial_conditions(mesh.get_point(i));
    weights.resize(mesh.get_number_of_points(), std::vector<double>((N+1)*(N+2)/2));
    //triangles_of_points.resize(mesh.get_number_of_points());
}



double gcmethod_2d::exact_solution( vector2d p, double t )
{
    vector2d r = p - t * vector2d(lambda_x, lambda_y);
    while (!mesh.is_inside(r))
        mesh.make_inside_vector(r);

    return initial_conditions(r);
}

double gcmethod_2d::L_inf()
{
    double max_err = fabs(values0.at(0) - exact_solution(mesh.get_point(0), tau*number_of_steps));
    for ( int i = 1; i < mesh.get_number_of_points(); i++ )
    {
        if ( max_err < fabs(values0.at(i) - exact_solution(mesh.get_point(i), tau*number_of_steps)) )
            max_err = fabs(values0.at(i) - exact_solution(mesh.get_point(i), tau*number_of_steps));
    }
    return max_err;
}

double gcmethod_2d::L(int n)
{
    double sum = 0, temp;
    for ( int i = 0; i < mesh.get_number_of_points(); i++ )
    {
        temp = pow(fabs(values0.at(i) - exact_solution(mesh.get_point(i), tau*number_of_steps)), n);
        if ( mesh.get_is_structured() )
        {
            temp *= mesh.get_h() * mesh.get_h();
            if (mesh.get_point(i).x < -mesh.get_size_x()/2 + eps_xy || mesh.get_point(i).x > mesh.get_size_x()/2 - eps_xy)
                temp /= 2;
            if (mesh.get_point(i).y < -mesh.get_size_y()/2 + eps_xy || mesh.get_point(i).y > mesh.get_size_y()/2 - eps_xy)
                temp /= 2;
            sum +=  temp;
        }
        else
        {
            sum += mesh.voronoi_areas.at(i) * temp;
        }
    }
    return pow(sum, 1.0/n) / mesh.get_size_x() / mesh.get_size_y();
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
    lambda_x = pt.get<REAL>("Equation.lambda_x");
    lambda_y = pt.get<REAL>("Equation.lambda_y");
    tau = pt.get<REAL>("Method.tau");
    number_of_steps = pt.get<int>("Method.number_of_steps");
    N = pt.get<int>("Method.N");
    std::string is_monotonic_str = pt.get<std::string>("Method.is_monotonic");
    is_monotonic = (is_monotonic_str == "true" || is_monotonic_str == "TRUE" || is_monotonic_str == "True");
}




