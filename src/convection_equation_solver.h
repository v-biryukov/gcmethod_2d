///////////////////////////////////////////////////////////
//////
//////	file: convection_equation_solver.h
//////  class for solving convection equation in 2d using
//////  grid-characteristic method.
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////


#pragma once

#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iomanip>
#include "vector2d.h"
#include "mesh_2d/mesh_2d.h"

struct convection_solver_gc_settings
{
    double lambda;
    vector2d direction;
    double tau;
    int number_of_steps;
    int N;
	bool is_monotonic;
    bool use_precalculations;
    bool save_only_main_points;
    int saving_frequency;
};


class gcmethod_2d
{
	std::string path;
	double lambda;
	vector2d direction;
	double tau;
	int number_of_steps;
	int cur_step;
    bool use_precalculations;
    bool save_only_main_points;
    int saving_frequency;

	double eps_xy = 1e-10;
	double eps_z = 1e-10;
	mesh_2d * mesh;

	std::vector<std::vector<double> > weights;
	std::vector<int> elements_of_points;

public:
    int N;
	bool is_monotonic;

    std::vector<double> values0;
	std::vector<double> values1;

    gcmethod_2d();
	gcmethod_2d(std::string path, mesh_2d * mesh_t);
	gcmethod_2d(struct convection_solver_gc_settings c, mesh_2d * mesh_t);
	void init();
	void calculate();
	void save_to_vtk(std::string name);
	void analyze();
    double L_inf();
    double L(int n);
    void set_number_of_steps(int n);
private:
	vector2d get_additional_point(int n, int k);
	void read_from_file(std::string path);
	double initial_conditions(double x, double y)
        {return 2*exp(-(x*x+y*y)/5);};//return 2*exp(-(x*x+y*y)/5);};//{if (x*x + y*y < 5) return 2.5; else return 0.0;};//return 2.0/(1 + x*x + y*y);};
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
    void fast_step();
    void calculate_weights(vector2d step);
    void calculate_weights_of_point(int pn, vector2d r, int tn);
    bool min_max_check(double z, int tn);
    double exact_solution(vector2d p, double t);


	void step();
};




