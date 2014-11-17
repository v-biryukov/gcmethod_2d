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

class mesh_2d
{
	// path to ini file
	std::string path;

	// sizes of the mesh
	double size_x, size_y;

	// area of the triangles < h*h/2
	double h, max_ang;

	// is mesh structured or unstructured
    bool is_structured;

	// number of divisions of horisontal and vertical borders
    int nx, ny;

    double eps = 1e-10;

	// step of the approximation
	// besides 3 main points of a triangle
	// (N+1)*(N+2)/2 - 3 additional points are creating for each element
    int N;

	// number of only main points (excluding additional)
	int number_of_main_points;

    // maximum triangle altitude
    double max_altitude;

public:

	// points of the mesh (including additional points)
	std::vector<vector2d> points;

	// elements of the mesh: (N+1)(N+2)/2 point numbers for each element
	std::vector<std::vector<int> > elements;

	// point neighbors of the point
	std::vector<std::vector<int> > neighbors;

	// triangle neighbors of the point
	std::vector<std::vector<int> > triangles;

    // maximum altitude among all elements
    double get_max_altitude();

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
    double get_size_x() {return size_x;};
	double get_size_y() {return size_y;};
    double get_h() {return h;};
	bool get_is_structured() {return is_structured;};
	vector2d get_point(int n);
	int get_opposite_point_num(int n);
	int get_number_of_triangles();

	//
	void change_h_and_refine(double ht);

	// get point/number of k's point of n's element
	int get_triangle_point_num(int n, int k);
	vector2d get_triangle_point(int n, int k);

	// check if p is inside mesh
	bool is_inside(vector2d p);
	// make p to be inside of the mesh (assuming mesh is periodical)
	void make_inside_vector(vector2d & p);
	// check if p is inside n's triangle
	bool is_inside(vector2d p, int n);


private:

	// write from triangulateio to points and elements (neighgours, triangles)
	void save_to_class_data(triangulateio * mesh);

	void read_from_file(std::string path);

	void find_neighbors(triangulateio * mesh);
	void find_triangles(triangulateio * mesh);
	void find_voronoi_areas();
    void find_max_altitude();

	//check if point n is corner point
	bool is_corner(int n);

	void create_structured_mesh();
	void create_unstructured_mesh();

    void init_triangulateio(triangulateio * in);
    void free_triangulateio(triangulateio * in);
};








