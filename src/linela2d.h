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

    void operator+=(point_data & pd)
    {
        vx += pd.vx;
        vy += pd.vy;
        sxx += pd.sxx;
        sxy += pd.sxy;
        syy += pd.syy;
    }

    point_data operator+(point_data & pd)
    {
        point_data temp;
        temp.vx = vx + pd.vx;
        temp.vy = vy + pd.vy;
        temp.sxx = sxx + pd.sxx;
        temp.sxy = sxy + pd.sxy;
        temp.syy = syy + pd.syy;
        return temp;
    }
};

struct riemann_data
{
    double w[5];
    riemann_data()
    {
        for (int k = 0; k < 5; k++)
            w[k] = 0;
    }
};

struct elements_data
{
    int el[2];
};

enum border_type {CONTINUOUS, ABSORB, SYMMETRIC};

class linela2d
{
    // pointer to a mesh
    mesh_2d * mesh;

    // Data on previous time layer
    std::vector<point_data> data;
    // Data on new time layer
    std::vector<point_data> data_new;

    // Indices of characteristic elements
    std::vector<elements_data> eldata;
    std::vector<elements_data> eldata_X;
    std::vector<elements_data> eldata_Y;

    // Parameters
    std::vector<double> c1;
    std::vector<double> c2;
    std::vector<double> rho;

    //
    std::vector<vector2d> directions;

    // Order of approximation
    int N;
    //
    int saving_frequency;
    double tau, h;
    int number_of_steps;
    bool is_monotonic;
    //
    bool is_axes_random;
    double lambda_hint;

    border_type hor_border_type = ABSORB;
    border_type vert_border_type = ABSORB;


    inline double get_rho(int n) {return rho[n];}
    inline double get_c1(int n) {return c1[n];}
    inline double get_c2(int n) {return c2[n];}
    inline double c3(int n) {return (c1[n]  - 2.0 * c2[n] * c2[n] / c1[n]);}
    inline double l(int n) {return rho[n] * (c1[n] * c1[n] - 2.0 * c2[n] * c2[n]);}
    inline double m(int n) {return rho[n] * c2[n] * c2[n];}

    std::vector<double> get_lambda_X(int point_n);
    std::vector<double> get_lambda_Y(int point_n);

    std::vector<int> get_point_elements(int pn);

    point_data rotate(point_data & origin, int sign);
    riemann_data get_riemann_inv_X(int point_n, int element_n);
    riemann_data get_riemann_inv_Y(int point_n, int element_n);
    riemann_data get_riemann_inv_X(int point_n, int element_n1, int element_n2);
    riemann_data get_riemann_inv_Y(int point_n, int element_n1, int element_n2);
    void set_point_data_X(int pn, riemann_data & rd);
    void set_point_data_Y(int pn, riemann_data & rd);

    void compute_rdata_X();
    void compute_rdata_Y();

    void calculate_point_elements();

    void step_X();
    void step_Y();
    void step();


    double approximate(vector2d r, int tn, std::vector<riemann_data> & rdata, int k);
    void set_directions(double angle);

    void read_from_file(std::string path);
    void init();


    void set_parameters()
    {
        for (int i = 0; i < mesh->get_number_of_triangles(); i++)
        {
            double x = (mesh->points[mesh->elements[i][0]].x + mesh->points[mesh->elements[i][1]].x + mesh->points[mesh->elements[i][2]].x)/3.0;
            double y = (mesh->points[mesh->elements[i][0]].y + mesh->points[mesh->elements[i][1]].y + mesh->points[mesh->elements[i][2]].y)/3.0;
            vector2d v = vector2d(x, y);


            int submesh_num = mesh->find_submesh(v);
            if (submesh_num < 0)
            {
                std::cout << "Error in Set parameters" << std::endl;
                std::exit(1);
            }
            else
            {
                c1[i] = mesh->get_c1(submesh_num);
                c2[i] = mesh->get_c2(submesh_num);
                rho[i] = mesh->get_rho(submesh_num);
            }

        }
    }

    void set_initial()
    {
        for (int i = 0; i < mesh->get_number_of_points(); i++)
        {
            double x = mesh->points[i].x;
            double y = mesh->points[i].y;
            vector2d v = vector2d(x, y);
            vector2d a = vector2d(0, 0);
            double start = 10;
            double finish = 30;
            if (x < finish && x > start && !mesh->triangles[i].empty())
            {
                int tn = mesh->triangles[i][0];
                data[i].vx = 1.0/(rho[tn]*c3(tn)) * sin((x-start)/(finish-start)*M_PI);
                data[i].vy = 0;
                data[i].sxx = c1[tn]/c3(tn) * sin((x-start)/(finish-start)*M_PI);
                data[i].sxy = 0;
                data[i].syy = 1 * sin((x-start)/(finish-start)*M_PI);
            }
            else
            {
                data[i].vx = 0;
                data[i].vy = 0;
                data[i].sxx = 0;
                data[i].sxy = 0;
                data[i].syy = 0;
            }
        }
    }

    // /////////////////////////////

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
