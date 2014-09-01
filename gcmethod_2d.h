
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

class mesh_2d
{
	std::string path;
	double size_x, size_y;
	double h, max_ang;
	
	triangulateio mesh;
	

public:
	mesh_2d(std::string path);
	void create_mesh();

	double get_size_x();
	double get_size_y();
	int get_number_of_points();
	vector2d get_point(int n);
	int get_number_of_triangles();
	int get_triangle_point_num(int n, int k);
	vector2d get_triangle_point(int n, int k);
	bool is_inside(vector2d p);
	bool is_inside(vector2d p, int n);

	std::vector<std::vector<int> > neighbors;
	std::vector<std::vector<int> > triangles;
private:
	void read_from_file(std::string path);
	void find_neighbors();
	void find_triangles();
};

mesh_2d::mesh_2d(std::string path)
{
	read_from_file(path);
}

void mesh_2d::read_from_file(std::string path)
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
    size_x = pt.get<REAL>("Grid.size_x");
    size_y = pt.get<REAL>("Grid.size_y");
    h = pt.get<REAL>("Method.h");
    max_ang = pt.get<REAL>("Method.max_ang");
}

void mesh_2d::create_mesh()
{
	// Setting triangulateio in and out:
	triangulateio in;
    in.numberofpoints = 4;
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    in.pointmarkerlist = (int *) NULL;
    in.pointattributelist = (REAL *) NULL;

    in.pointlist[0] = -size_x/2; in.pointlist[1] = -size_y/2;
    in.pointlist[2] =  size_x/2; in.pointlist[3] = -size_y/2;
    in.pointlist[4] =  size_x/2; in.pointlist[5] =  size_y/2;
    in.pointlist[6] = -size_x/2; in.pointlist[7] =  size_y/2;


    in.numberofsegments = in.numberofpoints;
    in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
	in.segmentlist[0] = 0; in.segmentlist[1] = 1;
	in.segmentlist[2] = 1; in.segmentlist[3] = 2;
	in.segmentlist[4] = 2; in.segmentlist[5] = 3;
	in.segmentlist[6] = 3; in.segmentlist[7] = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;

    mesh.pointlist = (REAL *) NULL;
    mesh.pointattributelist = (REAL *) NULL;
    mesh.pointmarkerlist = (int *) NULL;
    mesh.trianglelist = (int *) NULL;
    mesh.triangleattributelist = (REAL *) NULL;
    mesh.neighborlist = (int *) NULL;
    mesh.segmentlist = (int *) NULL;
    mesh.segmentmarkerlist = (int *) NULL;
    mesh.edgelist = (int *) NULL;
    mesh.edgemarkerlist = (int *) NULL;

    std::stringstream ss;
    ss << "pzeqa" << h * h / 2;
    std::string str;
    ss >> str;
    char * s = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), s);
    s[str.size()] = '\0';
    
    triangulate(s, &in, &mesh, (struct triangulateio *) NULL);
    find_neighbors();
    find_triangles();
    std::cout << "Mesh has been generated successfully." << std::endl;
}

void mesh_2d::find_neighbors()
{
	neighbors.resize(mesh.numberofpoints);
	for ( int i = 0; i < mesh.numberofedges; i++)
	{
		neighbors.at(mesh.edgelist[2*i + 0]).push_back(mesh.edgelist[2*i + 1]);
		neighbors.at(mesh.edgelist[2*i + 1]).push_back(mesh.edgelist[2*i + 0]);
	}
}

void mesh_2d::find_triangles()
{
	triangles.resize(mesh.numberofpoints);
	for ( int i = 0; i < mesh.numberoftriangles; i++)
	{
		triangles.at(mesh.trianglelist[3*i + 0]).push_back(i);
		triangles.at(mesh.trianglelist[3*i + 1]).push_back(i);
		triangles.at(mesh.trianglelist[3*i + 2]).push_back(i);
	}
}

double mesh_2d::get_size_x()
{
	return size_x;
}
double mesh_2d::get_size_y()
{
	return size_y;
}

int mesh_2d::get_number_of_points()
{
	return mesh.numberofpoints;
}

vector2d mesh_2d::get_point(int n)
{
	return vector2d(mesh.pointlist[2*n], mesh.pointlist[2*n+1]);
}

int mesh_2d::get_number_of_triangles()
{
	return mesh.numberoftriangles;
}

int mesh_2d::get_triangle_point_num(int n, int k)
{
	return mesh.trianglelist[3*n + k];
}

vector2d mesh_2d::get_triangle_point(int n, int k)
{
	return get_point(get_triangle_point_num(n,k));
}

bool mesh_2d::is_inside(vector2d p)
{
	if (p.x < -size_x/2 || p.x > size_x/2 || p.y < -size_y/2 || p.y > size_y/2)
		return false;
	else return true;
}

bool mesh_2d::is_inside(vector2d p, int n)
{
	bool b1, b2, b3;
	vector2d v1 = get_triangle_point(n, 0);
	vector2d v2 = get_triangle_point(n, 1);
	vector2d v3 = get_triangle_point(n, 2);
	b1 = Vec(p-v2, v1-v2);
	b2 = Vec(p-v3, v2-v3);
	b3 = Vec(p-v1, v3-v1);
	return ((b1 == b2) && (b2 == b3));
}




class gcmethod_2d
{	
	std::string path;
	double lambda_x, lambda_y;
	double tau;
	int number_of_steps;
	int N;
	mesh_2d & mesh;
	std::vector<double> u1;
	std::vector<double> u2;
	std::vector<std::vector<double> > values;

public:
	gcmethod_2d(std::string path, mesh_2d & mesh_t);
	void calculate();
	void save_to_vtk(std::string name);
private:
	void init();
	vector2d get_value_point(int n, int k);
	void read_from_file(std::string path);
	double initial_conditions(double x, double y) {if (x*x + y*y < 2) return 5; else return 0;};
	double initial_conditions(vector2d v) {return initial_conditions(v.x, v.y);};
	void step_any(vector2d step, std::vector<double> & u0, std::vector<double> & u);
	double approximate(vector2d p, int tn);
	void approximate_linear(vector2d p, int p1, int p2, int p3, std::vector<double> & u0, std::vector<double> & u);
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
		vtk_file << mesh.get_point(i).x << " " << mesh.get_point(i).y << " "  << u1.at(i) << "\n";
	vtk_file << "\nPOLYGONS " << mesh.get_number_of_triangles() << " " << mesh.get_number_of_triangles()*4 << "\n";
	for (int i = 0; i < mesh.get_number_of_triangles(); i++)
		vtk_file << 3 << " " << mesh.get_triangle_point_num(i,0) << " " << mesh.get_triangle_point_num(i,1) << " " << mesh.get_triangle_point_num(i,2) << "\n";
	vtk_file << "\nPOINT_DATA " << mesh.get_number_of_points() << "\n" << "VECTORS vectors float\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
		vtk_file << u1.at(i) << " 0.0 0.0\n";
}

void gcmethod_2d::step_any(vector2d step, std::vector<double> & u0, std::vector<double> & u)
{
	vector2d p;
	for (int i = 0; i < mesh.get_number_of_points(); i++)
	{
		p = mesh.get_point(i) + step;
		if (!mesh.is_inside(p))
			u.at(i) = 0;
		else
		{
			for (auto j : mesh.triangles.at(i))
			{
				if (mesh.is_inside(p, j))
				{
					u.at(i) = approximate(p, j);
					break;
				}
			}
			//approximate_linear(p, i, mesh.neighbors.at(i).at(k), mesh.neighbors.at(i).at(n), u0, u);
		}
	}
}

double gcmethod_2d::approximate(vector2d r, int tn)
{
	vector2d ra = mesh.get_triangle_point(tn, 0);
	vector2d rb = mesh.get_triangle_point(tn, 1);
	vector2d rc = mesh.get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb);
	double sb = Vec(ra - rc, r - rc);
	double sc = Vec(rc - ra, r - ra);
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
	double v = 0;
	for (int pnum = 0; pnum < values.at(tn).size(); pnum++)
	{
		int inum = 0, jnum = pnum, knum = 0;
		while ( jnum > inum )
		{
			inum++;
			jnum -= inum;
		}
		inum = N - inum;
		jnum = N - jnum;
		knum = N - inum - jnum;
		double w = 1;
		for (int i = 0; i < inum; i++)
			w *= (sa - i*1.0/N)/(1 - i*1.0/N);
		for (int j = 0; j < jnum; j++)
			w *= (sb - j*1.0/N)/(1 - j*1.0/N);
		for (int k = 0; k < knum; k++)
			w *= (sc - k*1.0/N)/(1 - k*1.0/N);
		v += w * values.at(tn).at(pnum);
	}
	return v;
}

void gcmethod_2d::approximate_linear(vector2d p, int p1, int p2, int p3, std::vector<double> & u0, std::vector<double> & u)
{
	double x1 = mesh.get_point(p1).x; double y1 = mesh.get_point(p1).y; double z1 = u0.at(p1);
	double x2 = mesh.get_point(p2).x; double y2 = mesh.get_point(p2).y; double z2 = u0.at(p2);
	double x3 = mesh.get_point(p3).x; double y3 = mesh.get_point(p3).y; double z3 = u0.at(p3);

	double vx = (y3-y1)*(z3-z2) - (z3-z1)*(y3-y2);
	double vy = (z3-z1)*(x3-x2) - (x3-x1)*(z3-z2);
	double vz = (x3-x1)*(y3-y2) - (y3-y1)*(x3-x2);

	u.at(p1) = 1.0 * (x1*vx + y1*vy + z1*vz - vx*p.x - vy*p.y) / vz ;
	
	// For deguging
	/*
	if ( (u.at(p1) > z1+1e-6 && u.at(p1) > z2+1e-6 && u.at(p1) > z3+1e-6) || vz == 0 )
	{
		std::cout << "Possible Error !!!\n";
		if (vz == 0) std::cout << "vz==0\n";
		else std::cout << "u>z, probably t * lambda is too big.\n";
		std::cout << "xy " << x1 << " " << y1 << " + " << x2 << " " << y2 << " + " << x3 << " " << y3 << "\n";
		std::cout << "pi " <<  p1 << " " << p2 <<  " " << p3 << "\n";
		std::cout << "p " <<  p.x << " " << p.y <<  " " << 0 << "\n";
		std::cout << "neighbors: \n";
		for (int i = 0; i < mesh.neighbors.at(p1).size(); i++)
			std::cout <<  mesh.neighbors.at(p1).at(i) << ": "  << " x: " << mesh.get_point(mesh.neighbors.at(p1).at(i)).x <<
			" y: " << mesh.get_point(mesh.neighbors.at(p1).at(i)).y << " z: " << u0.at(mesh.neighbors.at(p1).at(i)) << "\n";
		for (int i = 0; i < mesh.neighbors.at(p1).size(); i++)
			std::cout << mesh.get_point(mesh.neighbors.at(p1).at(i)).x << " " << mesh.get_point(mesh.neighbors.at(p1).at(i)).y <<"\n";
		std::cout << "v " <<  vx << " " << vy <<  " " << vz << "\n";
		std::cout << "z " <<  z1 << " " << z2 <<  " " << z3 << " " << u.at(p1) << "\n";
		//std::exit(1);
	}
	*/
	
}

void gcmethod_2d::step()
{
	// time step for the equation u_t + lx*u_x + ly*u_y = 0; 
	step_any(vector2d(-tau/2 * lambda_x, 0), u1, u2);
	step_any(vector2d(0, -tau/2 * lambda_y), u2, u1);
}

void gcmethod_2d::calculate()
{
	init();
	for (int i = 0; i < number_of_steps; i++)
	{
		// saving every step
		save_to_vtk("out/out_" + std::to_string(i) + ".vtk");
		step();
	}
	save_to_vtk("out/out_" + std::to_string(number_of_steps) + ".vtk");
}

void gcmethod_2d::init()
{
	u1.resize(mesh.get_number_of_points());
	u2.resize(mesh.get_number_of_points());
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
		u1.at(i) = initial_conditions(mesh.get_point(i).x, mesh.get_point(i).y);
	values.resize(mesh.get_number_of_triangles(), std::vector<double>((N+1)*(N+2)/2));
	for ( int i = 0; i < mesh.get_number_of_triangles(); i++ )
		for ( int j = 0; j < (N+1)*(N+2)/2; j++ )
			values.at(i).at(j) = initial_conditions(get_value_point(i, j));
}

vector2d gcmethod_2d::get_value_point(int n, int k)
{
	int i = 0;
	while ( k > i )
	{
		i++;
		k -= i;
	}
	vector2d ivec = (mesh.get_triangle_point(i, 2) - mesh.get_triangle_point(i, 0))/N;
	vector2d kvec = (mesh.get_triangle_point(i, 2) - mesh.get_triangle_point(i, 1))/N;
	return mesh.get_triangle_point(i, 0) + i*ivec + k*kvec;
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
}




