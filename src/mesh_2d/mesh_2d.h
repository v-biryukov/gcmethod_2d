///////////////////////////////////////////////////////////
//////
//////	file: mesh_2d.h
//////  class for 2d meshes using triangle c library
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

#define VOID int
#define ANSI_DECLARATORS


#if defined (__cplusplus)
extern "C" {
#endif

#include "triangle.h"

#if defined (__cplusplus)
}
#endif


#include "../vector2d.h"
#include "../tensor2d.h"
#include <string>
#include <vector>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

struct segment
{
    int start, finish;
};

struct submesh
{
    double rho, c1, c2;
    std::vector<segment> b_segments;
    std::vector<double> points_x;
    std::vector<double> points_y;
};

struct contour
{
    std::vector<segment> b_segments;
    int type;
};

class mesh_2d
{
	// path to ini file
	std::string path;

	// sizes of the mesh
    double size_x0, size_y0;
    double size_x1, size_y1;

    // submeshes
    std::string submeshes_file;
    std::vector<submesh> submeshes;

    // contours
    std::string contours_file;
    std::vector<point2d> cont_points;
    std::vector<contour> contours;

    // is mesh structured or unstructured
    bool is_complex;

	// area of the triangles < h*h/2
	double h, max_ang;

	// is mesh structured or unstructured
    bool is_structured;

	// number of divisions of horisontal and vertical borders
    int nx, ny;

    double eps;

	// step of the approximation
	// besides 3 main points of a triangle
	// (N+1)*(N+2)/2 - 3 additional points are creating for each element
    int N;

	// number of only main points (excluding additional)
	int number_of_main_points;

    // maximum triangle altitude
    double min_altitude;

    double courant_time_step;

public:

    enum point_type {INNER, FREE, FORCE, ABSORB};

    std::vector<point_type> border_sequence;

	// points of the mesh (including additional points)
	std::vector<vector2d> points;
    std::vector<point_type> point_types;
    std::vector<vector2d> point_normals;

	// elements of the mesh: (N+1)(N+2)/2 point numbers for each element
	std::vector<std::vector<int> > elements;

	// point neighbors of the point
	std::vector<std::vector<int> > neighbors;

	// triangle neighbors of the point
	std::vector<std::vector<int> > triangles;

    // maximum altitude among all elements
    double get_min_altitude();

	// voronoi areas for each point
	// area is restricted by borders
	std::vector<double> voronoi_areas;

	// CONSTRUCTOR
	mesh_2d(std::string path);

	// creating mesh
	void create_mesh();
	// Getters/setters
	int get_number_of_points();
	int get_number_of_main_points();
    double get_x0() {return size_x0;};
    double get_y0() {return size_x0;};
    double get_size_x() {return size_x1 - size_x0;};
    double get_size_y() {return size_y1 - size_y0;};
    double get_h() {return h;};
	bool get_is_structured() {return is_structured;};
	vector2d get_point(int n);
	int get_opposite_point_num(int n);
    int get_number_of_elements();
    int get_number_of_contour_points();

    vector2d get_center_point(int element_num);

	//
	void change_h_and_refine(double ht);

	// get point/number of k's point of n's element
	int get_triangle_point_num(int n, int k);
	vector2d get_triangle_point(int n, int k);

	// check if p is inside mesh
	bool is_inside(vector2d p);
	// make p to be inside of the mesh (assuming mesh is periodical)
    void make_inside_continuous(vector2d & p);
    // make p to be inside of the mesh (assuming border of the mesh is symmetric)
    void make_inside_symmetric(vector2d & p);
	// check if p is inside n's triangle
	bool is_inside(vector2d p, int n);
    // check if p is inside whole area
    bool is_inside_contour(vector2d p);
    // check if p is inside whole area (with eps)
    bool is_inside_contour(vector2d p, int pn);
    // find number of submesh point is in
    int find_submesh(vector2d p);

    double get_rho(int submesh_num) { return submeshes[submesh_num].rho; }
    double get_c1(int submesh_num) { return submeshes[submesh_num].c1; }
    double get_c2(int submesh_num) { return submeshes[submesh_num].c2; }

    double get_min_sound_speed();
    double get_max_sound_speed();

    void find_courant_time_step();

    double get_courant_time_step();

private:

	// write from triangulateio to points and elements (neighgours, triangles)
	void save_to_class_data(triangulateio * mesh);

	void read_from_file(std::string path);

	void find_neighbors(triangulateio * mesh);
	void find_triangles(triangulateio * mesh);
	void find_voronoi_areas();
    void find_min_altitude();

	//check if point n is corner point
	bool is_corner(int n);

    void load_submeshes(std::string spath);
    void load_contours(std::string cpath);
    void create_complex_mesh();


	void create_structured_mesh();
	void create_unstructured_mesh();

    int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy);

    void init_triangulateio(triangulateio * in);
    void free_triangulateio(triangulateio * in);
};








