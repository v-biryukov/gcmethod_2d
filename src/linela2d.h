///////////////////////////////////////////////////////////
//////
//////	file: linela2d.h
//////  class for solving linear elastisity equation in 2d using
//////  grid-characteristic method.
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////


#include "mesh_2d/mesh_2d.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <vector>

struct point_data
{
    double vx;
    double vy;
    double sxx;
    double sxy;
    double syy;
};

struct riemann_data
{
    double w[5];
};

struct elements_data
{
    int el[2];
};

class linela2d
{
    // pointer to a mesh
    mesh_2d * mesh;

    std::vector<point_data> data;
    //std::vector<point_data> data_prev;
    std::vector<riemann_data> rdata;
    std::vector<elements_data> eldata_X;
    std::vector<elements_data> eldata_Y;

    std::vector<double> c1;
    std::vector<double> c2;
    std::vector<double> rho;
    std::vector<vector2d> directions;

    int N;
    int saving_frequency;
    double tau, h;
    int number_of_steps;
    bool is_monotonic;
    bool is_axes_random;
    double lambda_hint;

    inline double get_rho(int n) {return rho[n];}
    inline double get_c1(int n) {return c1[n];}
    inline double get_c2(int n) {return c2[n];}
    inline double c3(int n) {return (c1[n]  - 2.0 * c2[n] * c2[n] / c1[n]);}
    inline double l(int n) {return rho[n] * (c1[n] * c1[n] - 2.0 * c2[n] * c2[n]);}
    inline double m(int n) {return rho[n] * c2[n] * c2[n];}

    std::vector<double> get_lambda_X(int point_n);
    std::vector<double> get_lambda_Y(int point_n);

    std::vector<int> get_point_elements(int pn);


    riemann_data get_riemann_inv_X(int point_n);
    riemann_data get_riemann_inv_Y(int point_n);
    void set_point_data_X(int pn, riemann_data & rd);
    void set_point_data_Y(int pn, riemann_data & rd);

    void compute_rdata_X();
    void compute_rdata_Y();

    void calculate_point_elements();

    void step_X();
    void step_Y();
    void step();


    double approximate(vector2d r, int tn, int k);
    void set_directions(double angle);

    void read_from_file(std::string path);
    void init();


    void set_parameters()
    {
        for (int i = 0; i < mesh->get_number_of_triangles(); i++)
        {
            c1[i] = 2.0; //* (mesh->points[mesh->elements[i][0]].x + 10.0)/10.0;
            c2[i] = 1.0; //* (mesh->points[mesh->elements[i][0]].x + 10.0)/10.0;
            rho[i] = 1.0;//(mesh->points[mesh->elements[i][0]].x + 10.0)/10.0;
        }
    }

    void set_initial()
    {
        for (int i = 0; i < mesh->get_number_of_points(); i++)
        {
            double x = mesh->points[i].x;
            double y = mesh->points[i].y;
            if (x*x + y*y < 1.5)
            {
                data[i].vx = (5 - x) / (x*x + y*y + 1);
                data[i].vy = (5 - y) / (x*x + y*y + 1);
            }
            /*if (y > -1 && y < 1)
            {
                data[i].vx = 0;
                data[i].vy = 1;
            }*/
            else
            {
                data[i].vx = 0;
                data[i].vy = 0;
            }
            data[i].sxx = 0;
            data[i].sxy = 0;
            data[i].syy = 0;
        }
    }

public:
    linela2d (std::string path, mesh_2d * m)
    {
        mesh = m;
        read_from_file(path);
        init();
        set_parameters();
        set_initial();
    }
    void calculate();
    void save_to_vtk(std::string name);
};
