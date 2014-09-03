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
#include "vector2.h"
#include "mesh_2d.h"






class gcmethod_2d
{
	std::string path;
	double lambda_x, lambda_y;
	double tau;
	int number_of_steps;
	int N;
	bool is_monotonic;
	mesh_2d & mesh;
	std::vector<double> main_z0;
	std::vector<double> main_z1;
	std::vector<std::vector<double> > additional_z0;
	std::vector<std::vector<double> > additional_z1;

public:
	gcmethod_2d(std::string path, mesh_2d & mesh_t);
	void calculate();
	void save_to_vtk(std::string name);
private:
	void init();
	vector2d get_additional_point(int n, int k);
	void read_from_file(std::string path);
	double initial_conditions(double x, double y) {if (x*x + y*y < 5) return 2; else return 0;};
	double initial_conditions(vector2d v) {return initial_conditions(v.x, v.y);};
	void step_any(vector2d step);
	double approximate(vector2d p, int tn);
	double approximate_linear(vector2d p, int tn);
	double approximate_quadratic(vector2d p, int tn);
    double approximate_cubic(vector2d p, int tn);
    double approximate_quartic(vector2d r, int tn);
    double approximate_any(vector2d p, int tn);
    bool min_max_check(double z, int tn);
	void step();
};

gcmethod_2d::gcmethod_2d(std::string path, mesh_2d & mesh_t) : mesh(mesh_t)
{
	read_from_file(path);
}

void gcmethod_2d::save_to_vtk(std::string name)
{
	std::ofstream vtk_file(name.c_str());
	vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
	vtk_file << "DATASET POLYDATA\nPOINTS " << mesh.get_number_of_points() << " float\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
		vtk_file << mesh.get_point(i).x << " " << mesh.get_point(i).y << " "  << main_z0.at(i) << "\n";
	vtk_file << "\nPOLYGONS " << mesh.get_number_of_triangles() << " " << mesh.get_number_of_triangles()*4 << "\n";
	for (int i = 0; i < mesh.get_number_of_triangles(); i++)
		vtk_file << 3 << " " << mesh.get_triangle_point_num(i,0) << " " << mesh.get_triangle_point_num(i,1) << " " << mesh.get_triangle_point_num(i,2) << "\n";
	vtk_file << "\nPOINT_DATA " << mesh.get_number_of_points() << "\n" << "VECTORS vectors float\n";
	//std::cout << "save points\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{
		//std::cout << values0.at(mesh.triangles.at(i).at(0)).at(0) << "\n";
		vtk_file <<  main_z0.at(i) << " 0.0 0.0\n";
	}

}

void gcmethod_2d::step_any(vector2d step)
{
	vector2d p;
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{
		p = mesh.get_point(i) + step;
		if (!mesh.is_inside(p))
			main_z1.at(i) = 0;
		else
		{
			for (auto n : mesh.triangles.at(i))
			{
				if (mesh.is_inside(p, n))
				{
					main_z1.at(i) = approximate(p, n);
					break;
				}
			}
		}
	}
	for ( int i = 0; i < mesh.get_number_of_triangles(); i++ )
	{
		for ( int j = 0; j < (N+1)*(N+2)/2 - 3; j++ )
		{
			p = get_additional_point(i, j) + step;
			if (!mesh.is_inside(p))
				additional_z1.at(i).at(j) = 0;
			else
			{
				for ( int m = 0; m < 3; m++ )
					for ( int n : mesh.triangles.at(mesh.get_triangle_point_num(i, m)) )
						if (mesh.is_inside(p, n))
						{
							additional_z1.at(i).at(j) = approximate(p, n);
							break;
						}
			}
		}
	}
}

bool gcmethod_2d::min_max_check(double z, int tn)
{
    return   (  z > main_z0.at(mesh.get_triangle_point_num(tn,0))
             && z > main_z0.at(mesh.get_triangle_point_num(tn,1))
             && z > main_z0.at(mesh.get_triangle_point_num(tn,2)) )
           ||(  z < main_z0.at(mesh.get_triangle_point_num(tn,0))
             && z < main_z0.at(mesh.get_triangle_point_num(tn,1))
             && z < main_z0.at(mesh.get_triangle_point_num(tn,2)) );

}

double gcmethod_2d::approximate(vector2d p, int tn)
{
    double result;
	switch (N)
	{
		case 1: result = approximate_linear(p, tn); break;
		case 2: result = approximate_quadratic(p, tn); break;
		case 3: result = approximate_cubic(p, tn); break;
		case 4: result = approximate_quartic(p, tn); break;
		default: result = approximate_any(p, tn); break;
	}
	if (is_monotonic && min_max_check(result, tn))
        return approximate_linear(p, tn);
    else
        return result;
}


double gcmethod_2d::approximate_linear(vector2d p, int tn)
{
	double x1 = mesh.get_triangle_point(tn, 0).x; double y1 = mesh.get_triangle_point(tn, 0).y; double z1 = main_z0.at(mesh.get_triangle_point_num(tn, 0));
	double x2 = mesh.get_triangle_point(tn, 1).x; double y2 = mesh.get_triangle_point(tn, 1).y; double z2 = main_z0.at(mesh.get_triangle_point_num(tn, 1));
	double x3 = mesh.get_triangle_point(tn, 2).x; double y3 = mesh.get_triangle_point(tn, 2).y; double z3 = main_z0.at(mesh.get_triangle_point_num(tn, 2));

	double vx = (y3-y1)*(z3-z2) - (z3-z1)*(y3-y2);
	double vy = (z3-z1)*(x3-x2) - (x3-x1)*(z3-z2);
	double vz = (x3-x1)*(y3-y2) - (y3-y1)*(x3-x2);

	return 1.0 * (x1*vx + y1*vy + z1*vz - vx*p.x - vy*p.y) / vz ;
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
	v += sa*(2*sa-1) * main_z0.at(mesh.get_triangle_point_num(tn, 0));
	v += sb*(2*sb-1) * main_z0.at(mesh.get_triangle_point_num(tn, 1));
	v += sc*(2*sc-1) * main_z0.at(mesh.get_triangle_point_num(tn, 2));
	v += 4*sa*sb * additional_z0.at(tn).at(0);
	v += 4*sc*sa * additional_z0.at(tn).at(1);
	v += 4*sb*sc * additional_z0.at(tn).at(2);
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
	v += sa*(3*sa-1)*(3*sa-2)/2 * main_z0.at(mesh.get_triangle_point_num(tn, 0));
	v += sb*(3*sb-1)*(3*sb-2)/2 * main_z0.at(mesh.get_triangle_point_num(tn, 1));
	v += sc*(3*sc-1)*(3*sc-2)/2 * main_z0.at(mesh.get_triangle_point_num(tn, 2));
	v += 9*sa*(3*sa-1)*sb/2 * additional_z0.at(tn).at(0);
	v += 9*sa*(3*sa-1)*sc/2 * additional_z0.at(tn).at(1);
	v += 9*sb*(3*sb-1)*sa/2 * additional_z0.at(tn).at(2);
	v += 27*sa*sb*sc * additional_z0.at(tn).at(3);
	v += 9*sc*(3*sc-1)*sa/2 * additional_z0.at(tn).at(4);
    v += 9*sb*(3*sb-1)*sc/2 * additional_z0.at(tn).at(5);
	v += 9*sc*(3*sc-1)*sb/2 * additional_z0.at(tn).at(6);
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
	v += sa*(4*sa-1)*(2*sa-1)*(4*sa-3)/3 * main_z0.at(mesh.get_triangle_point_num(tn, 0));
	v += sb*(4*sb-1)*(2*sb-1)*(4*sb-3)/3 * main_z0.at(mesh.get_triangle_point_num(tn, 1));
	v += sc*(4*sc-1)*(2*sc-1)*(4*sc-3)/3 * main_z0.at(mesh.get_triangle_point_num(tn, 2));
	v += 16*sa*(4*sa-1)*(2*sa-1)*sb/3 * additional_z0.at(tn).at(0);
	v += 16*sa*(4*sa-1)*(2*sa-1)*sc/3 * additional_z0.at(tn).at(1);

	v += 4*sa*(4*sa-1)*sb*(4*sb-1) * additional_z0.at(tn).at(2);
	v += 32*sa*(4*sa-1)*sb*sc * additional_z0.at(tn).at(3);
	v += 4*sa*(4*sa-1)*sc*(4*sc-1) * additional_z0.at(tn).at(4);
	v += 16*sb*(4*sb-1)*(2*sb-1)*sa/3 * additional_z0.at(tn).at(5);
    v += 32*sb*(4*sb-1)*sc*sa * additional_z0.at(tn).at(6);
    v += 32*sc*(4*sc-1)*sa*sb * additional_z0.at(tn).at(7);
    v += 16*sc*(4*sc-1)*(2*sc-1)*sa/3 * additional_z0.at(tn).at(8);
    v += 16*sb*(4*sb-1)*(2*sb-1)*sc/3 * additional_z0.at(tn).at(9);
    v += 4*sb*(4*sb-1)*sc*(4*sc-1) * additional_z0.at(tn).at(10);
    v += 16*sc*(4*sc-1)*(2*sc-1)*sb/3 * additional_z0.at(tn).at(11);
	return v;
}

double gcmethod_2d::approximate_any(vector2d r, int tn)
// TODO
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
	for (int pnum = 0; pnum < (N+1)*(N+2)/2; pnum++)
	{
		int inum = 0, jnum = 0, knum = pnum;
		while ( knum > inum )
		{
			inum++;
			knum -= inum;
		}
		inum = N - inum;
		jnum = N - inum - knum;
		double w = 1;
		//std::cout << inum*1.0/N;
		for (int i = 0; i < inum; i++)
            if (N==1)
                w *= sa;
            else
                w *= (sa - i*1.0/N)/(1 - i*1.0/N);
		for (int j = 0; j < jnum; j++)
            if (N==1)
                w *= sb;
            else
                w *= (sb - j*1.0/N)/(1 - j*1.0/N);
		for (int k = 0; k < knum; k++)
            if (N==1)
                w *= sc;
            else
                w *= (sc - k*1.0/N)/(1 - k*1.0/N);
		if (pnum == 0)
			v += w * main_z0.at(mesh.get_triangle_point_num(tn, 0));
		else if (pnum == (N+1)*(N+2)/2 - N - 1)
			v += w * main_z0.at(mesh.get_triangle_point_num(tn, 1));
		else if (pnum == (N+1)*(N+2)/2 - 1)
			v += w * main_z0.at(mesh.get_triangle_point_num(tn, 2));
		else
		{
            int k = pnum;
			k--;
			if ( k >= (N+1)*(N+2)/2 - N - 1 ) k--;
			v += w * additional_z0.at(tn).at(k);
		}
	}
	return v;
}

void gcmethod_2d::step()
{
	// time step for the equation u_t + lx*u_x + ly*u_y = 0;
	step_any(vector2d(-tau/2 * lambda_x, 0));
	main_z1.swap(main_z0);
	additional_z1.swap(additional_z0);
	step_any(vector2d(0, -tau/2 * lambda_y));
	main_z1.swap(main_z0);
	additional_z1.swap(additional_z0);
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
	main_z0.resize(mesh.get_number_of_points());
	main_z1.resize(mesh.get_number_of_points());
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
		main_z0.at(i) = initial_conditions(mesh.get_point(i));
    if (N>1)
    {
        additional_z1.resize(mesh.get_number_of_triangles(), std::vector<double>((N+1)*(N+2)/2 - 3));
        additional_z0.resize(mesh.get_number_of_triangles(), std::vector<double>((N+1)*(N+2)/2 - 3));
        for ( int i = 0; i < mesh.get_number_of_triangles(); i++ )
            for ( int j = 0; j < (N+1)*(N+2)/2 - 3; j++ )
            {
                additional_z0.at(i).at(j) = initial_conditions(get_additional_point(i, j));
                //std::cout << get_additional_point(i, j).x << " " << get_value_point(i, j).y << " " << std::endl;
            }
    }
}

vector2d gcmethod_2d::get_additional_point(int n, int k)
{
	k++;
	if (k >= (N+1)*(N+2)/2 - N - 1) k++;
	int i = 0;
	while ( k > i )
	{
		i++;
		k -= i;
	}
	vector2d ivec = (mesh.get_triangle_point(n, 1) - mesh.get_triangle_point(n, 0))/N;
	vector2d kvec = (mesh.get_triangle_point(n, 2) - mesh.get_triangle_point(n, 1))/N;
	return mesh.get_triangle_point(n, 0) + i*ivec + k*kvec;
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




