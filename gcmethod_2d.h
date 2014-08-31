
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

class gcmethod_2d
{	
	std::string path;
	double lambda_x, lambda_y;
	double size_x, size_y;
	double tau, h, max_ang;
	int number_of_steps;
	triangulateio mesh;
	std::vector<std::vector<int> > neighbors;
	std::vector<double> u1;
	std::vector<double> u2;

public:
	gcmethod_2d(std::string path);
	void create_mesh();
	void calculate();
	void save_to_vtk(std::string name);
private:
	void init();
	void read_from_file(std::string path);
	double initial_conditions(double x, double y) /*{return 1/(1 + x*x + y*y);};*/{if (x*x + y*y < 2) return 5; else return 0;};
	void find_neighbors(triangulateio & mesh);
	void step_any(vector2d step, std::vector<double> & u0, std::vector<double> & u);
	void approximate_linear(vector2d p, int p1, int p2, int p3, std::vector<double> & u0, std::vector<double> & u);
	void step();
};

gcmethod_2d::gcmethod_2d(std::string path)
{
	read_from_file(path);
}

void gcmethod_2d::save_to_vtk(std::string name)
{
	std::ofstream vtk_file(name.c_str());
	vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
	vtk_file << "DATASET POLYDATA\nPOINTS " << mesh.numberofpoints << " float\n";
	for ( int i = 0; i < mesh.numberofpoints; i++ )
		vtk_file << mesh.pointlist[2*i] << " " << mesh.pointlist[2*i+1] << " 0.0" << "\n";
	vtk_file << "\nPOLYGONS " << mesh.numberoftriangles << " " << mesh.numberoftriangles*4 << "\n";
	for (int i = 0; i < mesh.numberoftriangles; i++)
		vtk_file << 3 << " " << int(mesh.trianglelist[3*i]) << " " << int(mesh.trianglelist[3*i+1]) << " " << int(mesh.trianglelist[3*i+2]) << "\n";
	vtk_file << "\nPOINT_DATA " << mesh.numberofpoints << "\n" << "VECTORS vectors float\n";
	for ( int i = 0; i < mesh.numberofpoints; i++ )
		vtk_file << u1.at(i) << " 0.0 0.0\n";
}

void gcmethod_2d::step_any(vector2d step, std::vector<double> & u0, std::vector<double> & u)
{
	vector2d p0, p, q, r;
	for (int i = 0; i < mesh.numberofpoints; i++)
	{
		p0 = vector2d(mesh.pointlist[2*i + 0], mesh.pointlist[2*i + 1]);
		p = step;
		if ((p0 + p).x < -size_x/2 || (p0 + p).x > size_x/2 || (p0 + p).y < -size_y/2 || (p0 + p).y > size_y/2)
			u.at(i) = 0;
		else
		{
			double max_dot = -size_x - size_y;
			int k = 0;
			int n = 0;
			for ( int j = 0; j < neighbors.at(i).size(); j++ )
			{
				q = vector2d(mesh.pointlist[2*neighbors.at(i).at(j)]   - mesh.pointlist[2*i],
							 mesh.pointlist[2*neighbors.at(i).at(j)+1] - mesh.pointlist[2*i+1]);
				if (max_dot <= p*q/Magnitude(q)) 
				{
					max_dot = p*q/Magnitude(q);
					k = j;
				}
			}

			max_dot = -size_x - size_y;
			for ( int j = 0; j < neighbors.at(i).size(); j++ )
			{
				q = vector2d(mesh.pointlist[2*neighbors.at(i).at(j)]   - mesh.pointlist[2*i],
							 mesh.pointlist[2*neighbors.at(i).at(j)+1] - mesh.pointlist[2*i+1]);
				r = vector2d(mesh.pointlist[2*neighbors.at(i).at(k)]   - mesh.pointlist[2*i],
							 mesh.pointlist[2*neighbors.at(i).at(k)+1] - mesh.pointlist[2*i+1]);
				if (j != k && Vec(p, q) * Vec(p, r) <= 0)
				{
					if (max_dot <= p*q/Magnitude(q)) 
					{
						max_dot = p*q/Magnitude(q);
						n = j;
					}
				}
			}
			approximate_linear(p0 + p, i, neighbors.at(i).at(k), neighbors.at(i).at(n), u0, u);
		}
	}
}

void gcmethod_2d::approximate_linear(vector2d p, int p1, int p2, int p3, std::vector<double> & u0, std::vector<double> & u)
{
	double x1 = mesh.pointlist[2*p1 + 0]; double y1 = mesh.pointlist[2*p1 + 1]; double z1 = u0.at(p1);
	double x2 = mesh.pointlist[2*p2 + 0]; double y2 = mesh.pointlist[2*p2 + 1]; double z2 = u0.at(p2);
	double x3 = mesh.pointlist[2*p3 + 0]; double y3 = mesh.pointlist[2*p3 + 1]; double z3 = u0.at(p3);

	double vx = (y3-y1)*(z3-z2) - (z3-z1)*(y3-y2);
	double vy = (z3-z1)*(x3-x2) - (x3-x1)*(z3-z2);
	double vz = (x3-x1)*(y3-y2) - (y3-y1)*(x3-x2);

	u.at(p1) = 1.0 * (x1*vx + y1*vy + z1*vz - vx*p.x - vy*p.y) / vz ;
	/*
	// For deguging
	if ( (u.at(p1) > z1+1e-6 && u.at(p1) > z2+1e-6 && u.at(p1) > z3+1e-6) || vz == 0 )
	{
		std::cout << "Possible Error !!!\n";
		if (vz == 0) std::cout << "vz==0\n";
		else std::cout << "u>z, probably t * lambda is too big.\n";
		std::cout << "xy " << x1 << " " << y1 << " + " << x2 << " " << y2 << " + " << x3 << " " << y3 << "\n";
		std::cout << "pi " <<  p1 << " " << p2 <<  " " << p3 << "\n";
		std::cout << "p " <<  p.x << " " << p.y <<  " " << 0 << "\n";
		std::cout << "neighbors: \n";
		for (int i = 0; i < neighbors.at(p1).size(); i++)
			std::cout <<  neighbors.at(p1).at(i) << ": "  << " x: " << mesh.pointlist[2*neighbors.at(p1).at(i)] <<
			" y: " << mesh.pointlist[2*neighbors.at(p1).at(i)+1] << " z: " << u0.at(neighbors.at(p1).at(i)) << "\n";
		for (int i = 0; i < neighbors.at(p1).size(); i++)
			std::cout << mesh.pointlist[2*neighbors.at(p1).at(i)] << " " << mesh.pointlist[2*neighbors.at(p1).at(i)+1] <<"\n";
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
	//u1.swap(u2);
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
	u1.resize(mesh.numberofpoints);
	u2.resize(mesh.numberofpoints);
	for ( int i = 0; i < mesh.numberofpoints; i++ )
		u1.at(i) = initial_conditions(mesh.pointlist[2*i], mesh.pointlist[2*i+1]);
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
    size_x = pt.get<REAL>("Grid.size_x");
    size_y = pt.get<REAL>("Grid.size_y");
    lambda_x = pt.get<REAL>("Equation.lambda_x");
    lambda_y = pt.get<REAL>("Equation.lambda_y");
    tau = pt.get<REAL>("Method.tau");
    h = pt.get<REAL>("Method.h");
    max_ang = pt.get<REAL>("Method.max_ang");
    number_of_steps = pt.get<int>("Method.number_of_steps");
}

void gcmethod_2d::create_mesh()
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
    find_neighbors(mesh);
}

void gcmethod_2d::find_neighbors(triangulateio & mesh)
{
	neighbors.resize(mesh.numberofpoints);
	for ( int i = 0; i < mesh.numberofedges; i++)
	{
		neighbors.at(mesh.edgelist[2*i + 0]).push_back(mesh.edgelist[2*i + 1]);
		neighbors.at(mesh.edgelist[2*i + 1]).push_back(mesh.edgelist[2*i + 0]);
	}
	/* for debugging
	for ( int i = 0; i < mesh.numberofpoints; i++)
	{
		for ( int j = 0; j < neighbors.at(i).size(); j++)
			std::cout << neighbors.at(i).at(j) << " ";
		std::cout << "\n";
	}
	/* sort neighbors by angle
	for ( int i = 0; i < mesh.numberofpoints; i++)
	{
		int temp;
		for (int j = 0; j < neighbors.at(i).size() - 1; j++)
			for (int k = 0; k < neighbors.at(i).size() - j - 1; k++)
				if (atan2(mesh.pointlist[2*neighbors.at(i).at(k) + 1], mesh.pointlist[2*neighbors.at(i).at(k) + 0]) > 
					atan2(mesh.pointlist[2*neighbors.at(i).at(k+1) + 1], mesh.pointlist[2*neighbors.at(i).at(k+1) + 0]))
				{
					temp = neighbors.at(i).at(k);
					neighbors.at(i).at(k) = neighbors.at(i).at(k+1);
					neighbors.at(i).at(k+1) = neighbors.at(i).at(k);
				}
	}*/
}


