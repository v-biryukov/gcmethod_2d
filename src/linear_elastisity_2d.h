///////////////////////////////////////////////////////////
//////
//////	file: linear_elastisity_2d.h
//////  class for solving linear elastisity equation in 2d using
//////  grid-characteristic method.
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////

#pragma once

#include "vector2d.h"
#include "tensor2d.h"
#include "mesh_2d/mesh_2d.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

class linear_elastisity_2d
{
	// pointer to a mesh
	mesh_2d * mesh;

	// values
	std::vector<vector2d> values_v0;
	std::vector<tensor2d> values_T0;
	std::vector<vector2d> values_v1;
	std::vector<tensor2d> values_T1;

    double c1, c2, rho;

    std::vector<vector2d> xis;
    std::vector<std::vector<vector2d> > ns;
    std::vector<double> ls;
    std::vector<double> Lambda;
    std::vector<std::vector<tensor2d> > N0;
    std::vector<std::vector<tensor2d> > N1;

    int cur_step;
    double eps_z = 1e-10;


    int N;
	bool is_monotonic;
    bool use_precalculations;
    bool save_only_main_points;
    int saving_frequency;
    double tau, h;
	int number_of_steps;

    void init();
    void get_w(int j, std::vector<double> & w, vector2d v, tensor2d T);
    void get_vT (int j, std::vector<double> w, vector2d & v, tensor2d & T);
    void set_initial_conditions();
    void read_from_file(std::string path);
    void save_to_vtk(std::string name);

    double c3();
    void set_P_wave(vector2d & v, tensor2d & T, double x, double y, vector2d c, vector2d dir, double w, double mag);
    void set_S_wave(vector2d & v, tensor2d & T, double x, double y, vector2d c, vector2d dir, double w, double mag);
    void initial_conditions(vector2d & v, tensor2d & T, double x, double y);

    void initial_conditions(int j, std::vector<double> & w, double x, double y);
    double initial_conditions(int j, int k, double x, double y);
    void step();
    double approximate(std::vector<double> & w_prev, vector2d r, int tn);
    bool min_max_check(double z, std::vector<double> w_prev);
    double get_w(int j, int k, vector2d v, tensor2d T);

public:
    void calculate();
    linear_elastisity_2d (std::string path, mesh_2d * m)
    {
        mesh = m;
        read_from_file(path);
        init();
    }
};




