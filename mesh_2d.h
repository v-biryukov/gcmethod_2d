///////////////////////////////////////////////////////////
//////
//////  class for 2d meshes using triangle c library
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////

#pragma once

class mesh_2d
{
	std::string path;
	double size_x, size_y;
	double h, max_ang;
    bool is_structured;
    int nx, ny;
	triangulateio mesh;


public:
	mesh_2d(std::string path);
	void create_mesh();
	void create_structured_mesh();

	double get_size_x();
	double get_size_y();
	//int get_nx(){return static_cast<int>(size_x/h);};
	//int get_ny(){return static_cast<int>(size_y/h);};
	int get_number_of_points();
	vector2d get_point(int n);
	int get_opposite_point_num(int n);
	int get_number_of_triangles();
	int get_triangle_point_num(int n, int k);
	vector2d get_triangle_point(int n, int k);
	bool is_inside(vector2d p);
	bool is_inside(vector2d p, int n);
	void make_inside_vector(vector2d & p);

    double get_h() {return h;};
	bool get_is_structured() {return is_structured;};

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
    nx = static_cast<int>(size_x/h);
    ny = static_cast<int>(size_y/h);
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
    double step_x = size_x / nx;
    double step_y = size_y / ny;
	// Setting triangulateio in and out:
	triangulateio in;
    in.numberofpoints = 2*(nx + ny);
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    in.pointmarkerlist = (int *) NULL;
    in.pointattributelist = (REAL *) NULL;

    for ( int i = 0; i < nx; i++ )
    {
        in.pointlist[2*i + 0] = -size_x/2 + i * step_x;
        in.pointlist[2*i + 1] = -size_y/2;
        in.pointlist[2*(nx + ny) + 2*i + 0] = size_x/2 - i * step_x;
        in.pointlist[2*(nx + ny) + 2*i + 1] = size_y/2;
    }
    for ( int j = 0; j < ny; j++ )
    {
        in.pointlist[2*nx + 2*j + 0] = +size_x/2;
        in.pointlist[2*nx + 2*j + 1] = -size_y/2 + j * step_y;
        in.pointlist[2*(2*nx + ny) + 2*j + 0] = -size_x/2;
        in.pointlist[2*(2*nx + ny) + 2*j + 1] = +size_y/2 - j * step_y;
    }


    in.numberofsegments = in.numberofpoints;
    in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
    for ( int i = 0; i < in.numberofsegments-1; i++)
    {
        in.segmentlist[2*i + 0] = i;
        in.segmentlist[2*i + 1] = i+1;
    }
	in.segmentlist[2*(in.numberofsegments-1)] = in.numberofsegments-1;
	in.segmentlist[2*(in.numberofsegments-1) + 1] = 0;


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
    ss << "pzeYqa" << std::fixed << h * h / 2;
    std::string str;
    ss >> str;
    char * s = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), s);
    s[str.size()] = '\0';

    triangulate(s, &in, &mesh, (struct triangulateio *) NULL);
    find_neighbors();
    find_triangles();
    is_structured = false;
    std::cout << "Mesh has been generated successfully. " << std::endl;
}

void mesh_2d::create_structured_mesh()
{
    double step_x = size_x / nx;
    double step_y = size_y / ny;
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
    int pn = 0;
    for ( int i = 0; i < nx; i++ )
    {
        mesh.pointlist[pn++] = -size_x/2 + step_x * i;
        mesh.pointlist[pn++] = -size_y/2;
    }
    for ( int j = 0; j < ny; j++ )
    {
        mesh.pointlist[pn++] = +size_x/2;
        mesh.pointlist[pn++] = -size_y/2 + step_y * j;
    }
    for ( int i = 0; i < nx; i++ )
    {
        mesh.pointlist[pn++] = +size_x/2 - step_x * i;
        mesh.pointlist[pn++] = +size_y/2;
    }
    for ( int j = 0; j < ny; j++ )
    {
        mesh.pointlist[pn++] = -size_x/2;
        mesh.pointlist[pn++] = +size_y/2 - step_y * j;
    }
    for (int j = 1; j < ny; j++)
        for ( int i = 1; i < nx; i++ )
        {
            mesh.pointlist[pn++] = -size_x/2 + step_x * i;
            mesh.pointlist[pn++] = -size_y/2 + step_y * j;
        }
    int ei = 0, ti = 0;
    int p[4];
    for ( int j = 0; j < ny; j++ )
        for ( int i = 0; i < nx; i++ )
        {
            if ( j == 0 )
            {
                p[0] = i;
                p[1] = i+1;
                if ( i == nx-1 ) p[2] = nx+1; else p[2] = 2*(nx+ny) + i;
                p[3] = 2*(nx+ny) + i -1;
            }
            else if ( i == 0 )
            {
                p[0] = 2*(nx+ny) - j;
                p[1] = 2*(nx+ny) + (nx-1)*(j-1);
                if ( j == ny-1 ) p[2] = 2*(nx+ny) - j - 2; else p[2] = 2*(nx+ny) + (nx-1)*j;
                p[3] = 2*(nx+ny) - j - 1;
            }
            else if ( i == nx-1 )
            {
                p[0] = 2*(nx+ny) + nx - 2 + (nx-1)*(j-1);
                p[1] = nx + j;
                p[2] = nx + j + 1;
                if ( j == ny-1 ) p[3] = nx + j + 2; else p[3] = 2*(nx+ny) + nx - 2 + (nx-1)*j;
            }
            else if ( j == ny - 1 )
            {
                p[0] = 2*(nx+ny) + (nx-1)*(ny-2) + i - 1;
                p[1] = 2*(nx+ny) + (nx-1)*(ny-2) + i;
                p[2] = 2*nx + ny - i - 1;
                p[3] = 2*nx + ny - i;
            }
            else
            {
                p[0] = 2*(nx+ny) + (nx-1)*(j-1) + i - 1;
                p[1] = 2*(nx+ny) + (nx-1)*(j-1) + i;
                p[2] = 2*(nx+ny) + (nx-1)*(j) + i;
                p[3] = 2*(nx+ny) + (nx-1)*(j) + i - 1;
            }
            mesh.edgelist[ei++] =  p[0];
            mesh.edgelist[ei++] =  p[1];
            mesh.edgelist[ei++] =  p[0];
            mesh.edgelist[ei++] =  p[2];
            mesh.edgelist[ei++] =  p[0];
            mesh.edgelist[ei++] =  p[3];
            mesh.trianglelist[ti++] = p[0];
            mesh.trianglelist[ti++] = p[1];
            mesh.trianglelist[ti++] = p[2];
            mesh.trianglelist[ti++] = p[0];
            mesh.trianglelist[ti++] = p[2];
            mesh.trianglelist[ti++] = p[3];
            if ( i == nx - 1 )
            {
                mesh.edgelist[ei++] = p[1];
                mesh.edgelist[ei++] = p[2];
            }
            if ( j == ny - 1 )
            {
                mesh.edgelist[ei++] = p[3];
                mesh.edgelist[ei++] = p[2];
            }
        }
    find_neighbors();
    find_triangles();
    is_structured = true;
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
	for ( int i = 0; i < 2*(nx + ny); i++ )
        for ( auto x : triangles.at(get_opposite_point_num(i)) )
            if ( std::find(triangles.at(i).begin(), triangles.at(i).end(), x) == triangles.at(i).end() )
                triangles.at(i).push_back(x);

    std::vector<int> corners = {0, nx, nx+ny, 2*nx+ny};
    for (auto i : corners)
        for ( auto j : corners )
            for ( auto x : triangles.at(j) )
                if ( std::find(triangles.at(i).begin(), triangles.at(i).end(), x) == triangles.at(i).end() )
                    triangles.at(i).push_back(x);
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

int mesh_2d::get_opposite_point_num(int n)
{
    if ((n > 0) && (n < nx))
        return 2*nx + ny - n;
    else if ((n > nx) && (n < nx + ny))
        return 2*(nx+ny) - (n - nx);
    else if ((n > nx+ny) && (n < 2*nx + ny))
        return nx - (n-nx-ny);
    else if ((n > 2*nx+ny) && (n < 2*nx + 2*ny))
        return nx+ny - (n - 2*nx - ny);
    else if (n == 0 || n == nx || n == nx + ny || n == 2*nx + ny)
        return (n + nx + ny) % (2*nx + 2*ny);
    else
        std::cout << "Error in 'get_opposite_point_num'!\n";
}


