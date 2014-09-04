#pragma once

class mesh_2d
{
	std::string path;
	double size_x, size_y;
	double h, max_ang;

	triangulateio mesh;


public:
	mesh_2d(std::string path);
	void create_mesh();
	void create_structured_mesh();

	double get_size_x();
	double get_size_y();
	int get_number_of_points();
	vector2d get_point(int n);
	int get_number_of_triangles();
	int get_triangle_point_num(int n, int k);
	vector2d get_triangle_point(int n, int k);
	bool is_inside(vector2d p);
	bool is_inside(vector2d p, int n);
	void make_inside_vector(vector2d & p);

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
    ss << "pzeqa" << std::fixed << h * h / 2;
    std::string str;
    ss >> str;
    char * s = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), s);
    s[str.size()] = '\0';

    triangulate(s, &in, &mesh, (struct triangulateio *) NULL);
    find_neighbors();
    find_triangles();
    std::cout << "Mesh has been generated successfully. " << std::endl;
}

void mesh_2d::create_structured_mesh()
{
    int nx = static_cast<int>(size_x/h);
    int ny = static_cast<int>(size_y/h);
    mesh.numberofpoints = (nx+1)*(ny+1);
    mesh.pointlist = ((REAL *) malloc(mesh.numberofpoints * 2 * sizeof(REAL)));
    mesh.numberofpointattributes = 0;
    mesh.pointattributelist = (REAL *) NULL;
    mesh.pointmarkerlist = (int *) NULL;
    mesh.numberoftriangles = 2 * nx * ny;
    mesh.trianglelist = ((int *) malloc(mesh.numberoftriangles * 3 * sizeof(int)));
    mesh.numberoftriangleattributes = 0;
    mesh.triangleattributelist = (REAL *) NULL;
    mesh.neighborlist = (int *) NULL;
    mesh.segmentlist = (int *) NULL;
    mesh.segmentmarkerlist = (int *) NULL;
    mesh.numberofedges = 3 * nx * ny + nx + ny;
    mesh.edgelist = ((int *) malloc(mesh.numberofedges * 2 * sizeof(int)));
    mesh.edgemarkerlist = (int *) NULL;
    for (int j = 0; j < ny + 1; j++)
        for ( int i = 0; i < nx + 1; i++ )
        {
            mesh.pointlist[2*(j * (nx+1) + i) + 0] = -size_x/2 + h * i;
            mesh.pointlist[2*(j * (nx+1) + i) + 1] = -size_y/2 + h * j;
        }
    int ei = 0, ti = 0;
    for (int j = 0; j < ny; j++)
        for ( int i = 0; i < nx; i++ )
        {
            mesh.edgelist[ei++] =  j * (nx+1) + i;
            mesh.edgelist[ei++] =  j * (nx+1) + i + 1;
            mesh.edgelist[ei++] =  j * (nx+1) + i;
            mesh.edgelist[ei++] =  (j+1) * (nx+1) + i;
            mesh.edgelist[ei++] =  j * (nx+1) + i;
            mesh.edgelist[ei++] =  (j+1) * (nx+1) + i+1;
            mesh.trianglelist[ti++] = j * (nx+1) + i;
            mesh.trianglelist[ti++] = j * (nx+1) + i+1;
            mesh.trianglelist[ti++] = (j+1) * (nx+1) + i+1;
            mesh.trianglelist[ti++] = j * (nx+1) + i;
            mesh.trianglelist[ti++] = (j+1) * (nx+1) + i+1;
            mesh.trianglelist[ti++] = (j+1) * (nx+1) + i;
        }
    for (int j = 0; j < ny; j++)
    {
        mesh.edgelist[ei++] =  j * (nx+1) + nx;
        mesh.edgelist[ei++] =  (j+1) * (nx+1) + nx;
    }
    for ( int i = 0; i < nx; i++ )
    {
        mesh.edgelist[ei++] =  ny * (nx+1) + i;
        mesh.edgelist[ei++] =  ny * (nx+1) + i + 1;
    }
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

void mesh_2d::make_inside_vector(vector2d & p)
{
    if (p.x < -size_x/2) p.x += size_x;
    else if (p.x > size_x/2) p.x -= size_x;
    if (p.y < -size_y/2) p.y += size_y;
    else if (p.y > size_y/2) p.y -= size_y;
}

bool mesh_2d::is_inside(vector2d p, int n)
{
	bool b1, b2, b3;
	vector2d v1 = get_triangle_point(n, 0);
	vector2d v2 = get_triangle_point(n, 1);
	vector2d v3 = get_triangle_point(n, 2);
	b1 = Vec(p-v2, v1-v2) < 0.0;
	b2 = Vec(p-v3, v2-v3) < 0.0;
	b3 = Vec(p-v1, v3-v1) < 0.0;
	return ((b1 == b2) && (b2 == b3));
}


