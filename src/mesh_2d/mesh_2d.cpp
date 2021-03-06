///////////////////////////////////////////////////////////
//////
//////	file: mesh_2d.cpp
//////  class for 2d meshes using triangle c library
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////

#include "mesh_2d.h"
#include <fstream>
#include <sstream>
#include <iostream>

mesh_2d::mesh_2d(std::string path)
{
	read_from_file(path);
    nx = static_cast<int>((size_x1 - size_x0)/h);
    ny = static_cast<int>((size_y1 - size_y0)/h);
}

void mesh_2d::change_h_and_refine(double ht)
{
    h = ht;
    nx = static_cast<int>((size_x1 - size_x0)/h);
    ny = static_cast<int>((size_y1 - size_y0)/h);
	create_mesh();
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
    string is_complex_str = pt.get<std::string>("Mesh.is_complex");
    is_complex = ( is_complex_str == "true" || is_complex_str == "True" || is_complex_str == "TRUE" );

    std::stringstream ssx(pt.get<std::string>("Mesh.sizes_x"));
    ssx >> size_x0 >> size_x1;
    std::stringstream ssy(pt.get<std::string>("Mesh.sizes_y"));
    ssy >> size_y0 >> size_y1;

    h = pt.get<REAL>("Method.h");
    eps = 1e-6*h;
    max_ang = pt.get<REAL>("Method.max_ang");
    submeshes_file = pt.get<std::string>("Mesh.submeshes_file");
    contours_file = pt.get<std::string>("Mesh.contours_file");
    string is_structured_str = pt.get<std::string>("Mesh.is_structured");
    is_structured = ( is_structured_str == "true" || is_structured_str == "True" || is_structured_str == "TRUE" );

    N = pt.get<int>("Method.N");
    int number_of_border_types = pt.get<int>("Borders.number_of_border_types");
    border_sequence.push_back(INNER);
    for (int i = 1; i <= number_of_border_types; i++)
    {
        std::string current_type_str  = pt.get<std::string>("Borders.border_type_" + std::to_string(i));
        if (current_type_str == "ABSORB")
            border_sequence.push_back(ABSORB);
        else if (current_type_str == "FREE")
            border_sequence.push_back(FREE);
        else if (current_type_str == "FORCE")
            border_sequence.push_back(FORCE);
        else
        {
            std::cerr << "Error! There is no such border type " << current_type_str << "!" << std::endl;
            std::exit(1);
        }
    }
}
void mesh_2d::init_triangulateio(triangulateio * in)
{
    in->pointlist = (REAL*)(NULL);
    in->pointattributelist = (REAL*)(NULL);
    in->pointmarkerlist = (int*)(NULL);
    in->numberofpoints = 0;
    in->numberofpointattributes = 0;

    in->trianglelist = (int*)(NULL);
    in->triangleattributelist = (REAL*)(NULL);
    in->trianglearealist = (REAL*)(NULL);
    in->neighborlist = (int*)(NULL);
    in->numberoftriangles = 0;
    in->numberofcorners = 0;
    in->numberoftriangleattributes = 0;

    in->segmentlist = (int*)(NULL);
    in->segmentmarkerlist = (int*)(NULL);
    in->numberofsegments = 0;

    in->holelist = (REAL*)(NULL);
    in->numberofholes = 0;

    in->regionlist = (REAL*)(NULL);
    in->numberofregions = 0;

    in->edgelist = (int*)(NULL);
    in->edgemarkerlist = (int*)(NULL);
    in->normlist = (REAL*)(NULL);
    in->numberofedges = 0;
}

void mesh_2d::free_triangulateio(triangulateio * in)
{
    if (in->pointlist)             free(in->pointlist);
    if (in->pointattributelist)    free(in->pointattributelist);
    if (in->pointmarkerlist)       free(in->pointmarkerlist);
    if (in->trianglelist)          free(in->trianglelist);
    if (in->triangleattributelist) free(in->triangleattributelist);
    if (in->trianglearealist)      free(in->trianglearealist);
    if (in->neighborlist)          free(in->neighborlist);
    if (in->segmentlist)           free(in->segmentlist);
    if (in->segmentmarkerlist)     free(in->segmentmarkerlist);
    if (in->regionlist)            free(in->regionlist);
    if (in->edgelist)              free(in->edgelist);
    if (in->edgemarkerlist)        free(in->edgemarkerlist);
    if (in->normlist)              free(in->normlist);
}

int mesh_2d::pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

bool mesh_2d::is_inside_contour(vector2d p)
{
    return (find_submesh(p) >= 0);
}

bool mesh_2d::is_inside_contour(vector2d p, int pn)
{
    return (find_submesh(p-eps*point_normals[pn]) >= 0);
}

int mesh_2d::find_submesh(vector2d p)
{
    for (int i = 0; i < submeshes.size(); i++)
    {
        if (pnpoly(submeshes[i].b_segments.size(), submeshes[i].points_x.data(),  submeshes[i].points_y.data(),
                   p.x, p.y) == 1)
        {
            return i;
        }
    }
    return -1;
}

void mesh_2d::load_submeshes(std::string spath)
{
    std::ifstream is (spath);
    if (!is.is_open())
    {
        std::cout << "Unable to open file\n";
        std::exit(1);
    }
    int n_of_submeshes;
    is >> n_of_submeshes;
    submeshes.resize(n_of_submeshes);
    for (int i = 0; i < n_of_submeshes; i++)
    {
        int type;
        is >> type >> submeshes[i].rho >> submeshes[i].c1 >> submeshes[i].c2;
        int n_of_segments;
        is >> n_of_segments;
        submeshes[i].b_segments.resize(n_of_segments);
        submeshes[i].points_x.resize(n_of_segments);
        submeshes[i].points_y.resize(n_of_segments);
        for (int j = 0; j < n_of_segments; j++)
            is >> submeshes[i].b_segments[j].start >> submeshes[i].b_segments[j].finish;
    }
    int n_of_points;
    is >> n_of_points;
    std::vector<point2d> subm_points;
    subm_points.resize(n_of_points);
    for (int i = 0; i < n_of_points; i++)
    {
        is >> subm_points[i].x >> subm_points[i].y;
    }
    for (int i = 0; i < n_of_submeshes; i++)
    {
        for (int j = 0; j < submeshes[i].b_segments.size(); j++)
        {
            submeshes[i].points_x[j] = subm_points[submeshes[i].b_segments[j].start].x;
            submeshes[i].points_y[j] = subm_points[submeshes[i].b_segments[j].start].y;
        }
    }
}


void mesh_2d::load_contours(std::string cpath)
{
    std::ifstream is (cpath);
    if (!is.is_open())
    {
        std::cout << "Unable to open file\n";
        std::exit(1);
    }
    int n_of_contours;
    is >> n_of_contours;
    contours.resize(n_of_contours);
    for (int i = 0; i < n_of_contours; i++)
    {
        double c;
        int n_of_segments;
        int type_num;
        is >> type_num >> c >> c >> c >> n_of_segments;
        contours[i].type = (int)border_sequence[type_num];
        contours[i].b_segments.resize(n_of_segments);
        for (int j = 0; j < n_of_segments; j++)
            is >> contours[i].b_segments[j].start >> contours[i].b_segments[j].finish;
    }
    int n_of_points;
    is >> n_of_points;
    cont_points.resize(n_of_points);
    is >> cont_points[0].x >> cont_points[0].y;
    double max_x = cont_points[0].x;
    double max_y = cont_points[0].y;
    double min_x = cont_points[0].x;
    double min_y = cont_points[0].y;
    for (int i = 1; i < n_of_points; i++)
    {
        is >> cont_points[i].x >> cont_points[i].y;
        if (cont_points[i].x > max_x) max_x = cont_points[i].x;
        if (cont_points[i].x < min_x) min_x = cont_points[i].x;
        if (cont_points[i].y > max_y) max_y = cont_points[i].y;
        if (cont_points[i].y < min_y) min_y = cont_points[i].y;
    }
    size_x0 = min_x;
    size_x1 = max_x;
    size_y0 = min_y;
    size_y1 = max_y;
}

void mesh_2d::create_mesh()
{
    if (is_complex)
    {
        load_submeshes(submeshes_file);
        load_contours(contours_file);
        create_complex_mesh();
        find_min_altitude();
        find_courant_time_step();
    }
    else
    {
        if ( is_structured )
            create_structured_mesh();
        else
            create_unstructured_mesh();
        find_min_altitude();
    }
}

void mesh_2d::create_complex_mesh()
{
    triangulateio in, mesh;
    init_triangulateio(&in);
    init_triangulateio(&mesh);
    in.numberofpoints = cont_points.size();
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    in.pointmarkerlist = (int *) NULL;
    in.pointattributelist = (REAL *) NULL;
    for (int i = 0; i < in.numberofpoints; i++)
    {
        in.pointlist[2*i] = cont_points[i].x;
        in.pointlist[2*i+1] = cont_points[i].y;
    }
    in.numberofsegments = 0;
    for (int i = 0; i < contours.size(); i++)
    {
        in.numberofsegments += contours[i].b_segments.size();
    }
    in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
    int k = 0;
    for (int i = 0; i < contours.size(); i++)
    {
        for (int j = 0; j < contours[i].b_segments.size(); j++)
        {
            in.segmentlist[2*k] = contours[i].b_segments[j].start;
            in.segmentlist[2*k+1] = contours[i].b_segments[j].finish;
            k += 1;
        }
    }
    in.numberofholes = 0;
    in.numberofregions = 0;

    std::stringstream ss;
    ss << "pzeQYqa" << std::fixed << h * h / 2;
    std::string str;
    ss >> str;
    char * s = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), s);
    s[str.size()] = '\0';
    triangulate(s, &in, &mesh, (struct triangulateio *) NULL);


    // Setting point types and normals
    point_types.resize(mesh.numberofpoints, INNER);
    point_normals.resize(mesh.numberofpoints, vector2d(0.0, 0.0));
    for (int i = 0; i < contours.size(); i++)
    {
        for (int j = 0; j < contours[i].b_segments.size(); j++)
        {
            // Some edges can be trivial. TODO
            if (contours[i].b_segments[j].start == contours[i].b_segments[j].finish)
                continue;
            vector2d normal = (cont_points[contours[i].b_segments[j].finish]-cont_points[contours[i].b_segments[j].start]).rotate(M_PI/2);
            normal = normal.Normalize();
            vector2d middle = (cont_points[contours[i].b_segments[j].finish]+cont_points[contours[i].b_segments[j].start])/2.0;
            if (is_inside_contour(middle + normal*h*0.1))
                normal *= -1.0;
            if ((int)point_types[contours[i].b_segments[j].start] < contours[i].type)
            {
                point_types[contours[i].b_segments[j].start] = (point_type)contours[i].type;
                point_normals[contours[i].b_segments[j].start]  = normal;
            }
            if ((int)point_types[contours[i].b_segments[j].finish] < contours[i].type)
            {
                point_types[contours[i].b_segments[j].finish] = (point_type)contours[i].type;
                point_normals[contours[i].b_segments[j].finish] = normal;
            }

            //point_normals[contours[i].b_segments[j].start]  += normal;
            //point_normals[contours[i].b_segments[j].finish] += normal;

        }
    }

    for (int i = 0; i < point_normals.size(); i++)
    {
        if (Magnitude(point_normals[i]) > 0.1) // TODO
        {
            point_normals[i] = point_normals[i].Normalize();
        }
    }

    save_to_class_data(&mesh);
    free_triangulateio(&in);
    free_triangulateio(&mesh);
    std::cout << "Mesh has been generated successfully." << std::endl;
}

void mesh_2d::create_unstructured_mesh()
{
    bool is_line = false;
    double step_x = (size_x1 - size_x0) / nx;
    double step_y = (size_y1 - size_y0) / ny;
	// Setting triangulateio in and out:
	triangulateio in, mesh;
    init_triangulateio(&in);
    init_triangulateio(&mesh);
    in.numberofpoints = 2*(nx + ny);
    if (is_line)
        in.numberofpoints += ny - 1;
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    in.pointmarkerlist = (int *) NULL;
    in.pointattributelist = (REAL *) NULL;

    for ( int i = 0; i < nx; i++ )
    {
        in.pointlist[2*i + 0] = size_x0 + i * step_x;
        in.pointlist[2*i + 1] = size_y0;
        in.pointlist[2*(nx + ny) + 2*i + 0] = size_x1 - i * step_x;
        in.pointlist[2*(nx + ny) + 2*i + 1] = size_y1;
    }
    for ( int j = 0; j < ny; j++ )
    {
        in.pointlist[2*nx + 2*j + 0] = +size_x1;
        in.pointlist[2*nx + 2*j + 1] = size_y0 + j * step_y;
        in.pointlist[2*(2*nx + ny) + 2*j + 0] = size_x0;
        in.pointlist[2*(2*nx + ny) + 2*j + 1] = +size_y1 - j * step_y;
    }
    if (is_line)
    {
        for ( int j = 1; j < ny; j++ )
        {
            in.pointlist[2*(2*nx + 2*ny) + 2*(j-1) + 0] = 0.0;
            in.pointlist[2*(2*nx + 2*ny) + 2*(j-1) + 1] = size_y0 + j * step_y;
        }
    }


    in.numberofsegments = in.numberofpoints;
    if (is_line)
        in.numberofsegments += 1;
    in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
    for ( int i = 0; i < 2*(nx + ny)-1; i++)
    {
        in.segmentlist[2*i + 0] = i;
        in.segmentlist[2*i + 1] = i+1;
    }
    in.segmentlist[2*(2*(nx + ny)-1)] = 2*(nx + ny)-1;
    in.segmentlist[2*(2*(nx + ny)-1) + 1] = 0;

    if (is_line)
    {
        in.segmentlist[2*(2*(nx + ny))] = nx / 2;
        in.segmentlist[2*(2*(nx + ny)) + 1] = 2*(nx + ny);
        for ( int i = 0; i < ny-2; i++)
        {
            in.segmentlist[2*(2*(nx + ny)+i) + 0] = 2*(nx + ny) + i;
            in.segmentlist[2*(2*(nx + ny)+i) + 1] = 2*(nx + ny) + i+1;
        }
        in.segmentlist[2*(2*(nx + ny))] = 2*nx + ny - nx / 2;
        in.segmentlist[2*(2*(nx + ny)) + 1] = in.numberofpoints - 1;
    }

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
    ss << "pzeQYqa" << std::fixed << h * h / 2;
    std::string str;
    ss >> str;
    char * s = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), s);
    s[str.size()] = '\0';

    triangulate(s, &in, &mesh, (struct triangulateio *) NULL);
    save_to_class_data(&mesh);
    free_triangulateio(&in);
    free_triangulateio(&mesh);
    std::cout << "Mesh has been generated successfully." << std::endl;
}

void mesh_2d::save_to_class_data(triangulateio * mesh)
{
    // Main points
    number_of_main_points = mesh->numberofpoints;
    for ( int i = 0; i < mesh->numberofpoints; i++ )
        points.push_back(vector2d(mesh->pointlist[2*i+0], mesh->pointlist[2*i+1]));
    int points_per_el = (N+1)*(N+2)/2;
    elements.resize(mesh->numberoftriangles, std::vector<int>(points_per_el));
    for ( int i = 0; i < mesh->numberoftriangles; i++ )
    {
        elements.at(i).at(0) = mesh->trianglelist[3*i + 0];
        elements.at(i).at(1) = mesh->trianglelist[3*i + 1];
        elements.at(i).at(2) = mesh->trianglelist[3*i + 2];
    }
    find_neighbors(mesh);
    find_triangles(mesh);
    // Additional points:
    for ( int i = 0; i < mesh->numberofedges; i++ )
    {
        vector2d dp = (points.at(mesh->edgelist[2*i + 1]) - points.at(mesh->edgelist[2*i + 0]))/N;
        for ( int k = 1; k < N; k++ )
        {
            points.push_back(points.at(mesh->edgelist[2*i + 0]) + k * dp);
            point_types.push_back(std::min(point_types.at(mesh->edgelist[2*i + 1]), point_types.at(mesh->edgelist[2*i + 0])));
            vector2d normal;
            if (int(point_types[mesh->edgelist[2*i + 0]]) > int(point_types[mesh->edgelist[2*i + 1]]))
                normal = point_normals[mesh->edgelist[2*i + 1]];
            else
                normal = point_normals[mesh->edgelist[2*i + 0]];
            point_normals.push_back(normal);
            int ep0 = mesh->edgelist[2*i + 0];
            int ep1 = mesh->edgelist[2*i + 1];
            for (auto t : triangles.at(ep0))
                for (auto s : triangles.at(ep1))
                {
                    if (t == s)
                        {
                            std::vector<int> & ts = elements.at(t);
                            if (ts.at(0) == ep0 && ts.at(1) == ep1) ts.at(3 + k - 1) = points.size()-1;
                            else if (ts.at(1) == ep0 && ts.at(0) == ep1) ts.at(3 + N - k - 1) = points.size()-1;
                            else if (ts.at(1) == ep0 && ts.at(2) == ep1) ts.at(3 + N + k - 2) = points.size()-1;
                            else if (ts.at(2) == ep0 && ts.at(1) == ep1) ts.at(3 + 2*N - k - 2) = points.size()-1;
                            else if (ts.at(2) == ep0 && ts.at(0) == ep1) ts.at(3 + 2*N + k - 3) = points.size()-1;
                            else if (ts.at(0) == ep0 && ts.at(2) == ep1) ts.at(3 + 3*N - k - 3) = points.size()-1;
                        }
                }
            triangles.push_back(triangles.at(mesh->edgelist[2*i + 0]));
            triangles.back().insert(triangles.back().end(), triangles.at(mesh->edgelist[2*i + 1]).begin(), triangles.at(mesh->edgelist[2*i + 1]).end());
        }
    }
    for ( int i = 0; i < mesh->numberoftriangles; i++ )
    {
        for ( int k = 0; k < (N-1)*(N-2)/2; k++ )
        {
            int l = 0, m = k;
            while ( m > l )
            {
                l++;
                m -= l;
            }
            vector2d lvec = (points.at(mesh->trianglelist[3*i + 1]) - points.at(mesh->trianglelist[3*i + 0]))/N;
            vector2d mvec = (points.at(mesh->trianglelist[3*i + 2]) - points.at(mesh->trianglelist[3*i + 1]))/N;
            points.push_back(points.at(mesh->trianglelist[3*i + 0]) + (l+2)*lvec + (m+1)*mvec);
            point_types.push_back(INNER);
            elements.at(i).at(3 + 3*(N-1) + k) = points.size()-1;
            triangles.push_back(triangles.at(mesh->trianglelist[3*i + 0]));
            triangles.back().insert(triangles.back().end(), triangles.at(mesh->trianglelist[3*i + 1]).begin(), triangles.at(mesh->trianglelist[3*i + 1]).end());
            triangles.back().insert(triangles.back().end(), triangles.at(mesh->trianglelist[3*i + 2]).begin(), triangles.at(mesh->trianglelist[3*i + 2]).end());
        }
    }
    if (false)//!is_complex)
        find_voronoi_areas();
}

void mesh_2d::create_structured_mesh()
{
    double step_x = (size_x1 - size_x0) / nx;
    double step_y = (size_y1 - size_y0) / ny;
    triangulateio mesh;
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
        mesh.pointlist[pn++] = size_x0 + step_x * i;
        mesh.pointlist[pn++] = size_y0;
    }
    for ( int j = 0; j < ny; j++ )
    {
        mesh.pointlist[pn++] = +size_x1;
        mesh.pointlist[pn++] = size_y0 + step_y * j;
    }
    for ( int i = 0; i < nx; i++ )
    {
        mesh.pointlist[pn++] = +size_x1 - step_x * i;
        mesh.pointlist[pn++] = +size_y1;
    }
    for ( int j = 0; j < ny; j++ )
    {
        mesh.pointlist[pn++] = size_x0;
        mesh.pointlist[pn++] = +size_y1 - step_y * j;
    }
    for (int j = 1; j < ny; j++)
        for ( int i = 1; i < nx; i++ )
        {
            mesh.pointlist[pn++] = size_x0 + step_x * i;
            mesh.pointlist[pn++] = size_y0 + step_y * j;
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
    save_to_class_data(&mesh);
    free(mesh.pointlist);
    free(mesh.segmentlist);
    free(mesh.edgelist);
    free(mesh.trianglelist);
    std::cout << "Mesh has been generated successfully." << std::endl;
}

void mesh_2d::find_neighbors(triangulateio * mesh)
{
	neighbors.resize(mesh->numberofpoints);
	for ( int i = 0; i < mesh->numberofedges; i++)
	{
		neighbors.at(mesh->edgelist[2*i + 0]).push_back(mesh->edgelist[2*i + 1]);
		neighbors.at(mesh->edgelist[2*i + 1]).push_back(mesh->edgelist[2*i + 0]);
	}
}

void mesh_2d::find_triangles(triangulateio * mesh)
{
	triangles.resize(mesh->numberofpoints);
	for ( int i = 0; i < mesh->numberoftriangles; i++)
	{
		triangles.at(mesh->trianglelist[3*i + 0]).push_back(i);
		triangles.at(mesh->trianglelist[3*i + 1]).push_back(i);
		triangles.at(mesh->trianglelist[3*i + 2]).push_back(i);
	}
    if (!is_complex)
    {
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
}


void mesh_2d::find_min_altitude()
{
    min_altitude = size_x1 - size_x0 + size_y1 - size_y0;
    for (int i = 0; i < get_number_of_elements(); i++)
    {
        vector2d ra = points[elements[i][1]] - points[elements[i][0]];
        vector2d rb = points[elements[i][2]] - points[elements[i][1]];
        vector2d rc = points[elements[i][0]] - points[elements[i][2]];
        double ha = Magnitude(ra - Dot(ra, rb) * rb / SquaredMag(rb) );
        double hb = Magnitude(rb - Dot(rb, rc) * rc / SquaredMag(rc) );
        double hc = Magnitude(rc - Dot(rc, ra) * ra / SquaredMag(ra) );
        //std::cout << i << " : " << min_altitude << " : " << ha << " : " << hb << " : " << hc <<std::endl;
        if (ha < min_altitude) min_altitude = ha;
        if (hb < min_altitude) min_altitude = hb;
        if (hc < min_altitude) min_altitude = hc;
    }
}



double mesh_2d::get_min_altitude()
{
    return min_altitude;
}

double mesh_2d::get_min_sound_speed()
{
    double min_sound_speed = std::min(submeshes[0].c1, submeshes[0].c2);
    for (int i = 0; i < submeshes.size(); i++)
    {
        if (submeshes[i].c1 < min_sound_speed) min_sound_speed = submeshes[i].c1;
        if (submeshes[i].c2 < min_sound_speed) min_sound_speed = submeshes[i].c2;
    }
    return min_sound_speed;
}

double mesh_2d::get_max_sound_speed()
{
    double max_sound_speed = std::max(submeshes[0].c1, submeshes[0].c2);
    for (int i = 0; i < submeshes.size(); i++)
    {
        if (submeshes[i].c1 > max_sound_speed) max_sound_speed = submeshes[i].c1;
        if (submeshes[i].c2 > max_sound_speed) max_sound_speed = submeshes[i].c2;
    }
    return max_sound_speed;
}

void mesh_2d::find_courant_time_step()
{
    courant_time_step = 10.0 * h / get_min_sound_speed();
    for (int i = 0; i < get_number_of_elements(); i++)
    {
        vector2d element_center = (points[elements[i][0]]+points[elements[i][1]]+points[elements[i][2]])/3.0;
        int submesh_num = find_submesh(element_center);
        double c;
        if (submesh_num < 0)
            c = get_max_sound_speed();
        else
            c = std::max(submeshes.at(submesh_num).c1, submeshes.at(submesh_num).c2);//=c1

        vector2d ra = points[elements[i][1]] - points[elements[i][0]];
        vector2d rb = points[elements[i][2]] - points[elements[i][1]];
        vector2d rc = points[elements[i][0]] - points[elements[i][2]];
        double ha = Magnitude(ra - Dot(ra, rb) * rb / SquaredMag(rb) );
        double hb = Magnitude(rb - Dot(rb, rc) * rc / SquaredMag(rc) );
        double hc = Magnitude(rc - Dot(rc, ra) * ra / SquaredMag(ra) );

        double courant_a = ha / c / N;
        double courant_b = hb / c / N;
        double courant_c = hc / c / N;

        if (courant_a < courant_time_step) courant_time_step = courant_a;
        if (courant_b < courant_time_step) courant_time_step = courant_b;
        if (courant_c < courant_time_step) courant_time_step = courant_c;
    }
}

double mesh_2d::get_courant_time_step()
{
    return courant_time_step;
}

// Function that finds voronoi areas for each point and saves
// voronoi diagram in out/voronoi.vtk (only for simple rect mesh)
void mesh_2d::find_voronoi_areas()
{
    // Create file
    std::ofstream vtk_file("out/voronoi.vtk");
    vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
    std::vector<vector2d> voronp;
    std::vector<std::vector<int> > voront;
    ////////////////////////////////
    voronoi_areas.resize(points.size());
    // coefs - vector of lines - 1 for each neighbor ( +1 for each neighbor of the neighbors ) +1 for each border
    std::vector< std::vector<double> > coefs;
    // voronoi_local_points - points, that are constitute voronoi cell of the point i
    std::vector<vector2d> voronoi_local_points;
    for ( int i = 0; i < number_of_main_points; i++ )
    {
        // lines which correspond to the borders
        coefs.push_back({0, 1, size_y0});
        coefs.push_back({0, -1,size_y0});
        coefs.push_back({1, 0, size_x0});
        coefs.push_back({-1, 0, size_x0});
        // lines which correspond to the neighbors
        for (int j : neighbors.at(i))
        {
            vector2d p = get_point(i);
            vector2d q = get_point(j);
            std::vector<double> temp = {(q-p).x, (q-p).y, (p*p - q*q)/2};
            coefs.push_back(temp);
        }

        // lines which correspond to the neighbors of the neighbors
        for (int j : neighbors.at(i))
        {
            for (int k : neighbors.at(j))
            {
                if (k != i)
                {
                    vector2d p = get_point(i);
                    vector2d q = get_point(k);
                    std::vector<double> temp = {(q-p).x, (q-p).y, (p*p - q*q)/2};
                    bool is_in_coefs = false;
                    for (auto c : coefs)
                        if (fabs(c[0] - temp[0]) < eps && fabs(c[1] - temp[1]) < eps && fabs(c[2] - temp[2]) < eps)
                        {
                            is_in_coefs = true;
                            break;
                        }
                    if ( !is_in_coefs ) coefs.push_back(temp);
                }
            }
        }

        // In the next cycle we are looking for intersections of the lines
        // which are lying to the negative side of each line (or lie on line)
        for ( unsigned int m = 0; m < coefs.size(); m++ )
            for ( unsigned int n = m+1; n < coefs.size(); n++ )
            {
                double det = (coefs.at(m).at(1) * coefs.at(n).at(0) - coefs.at(m).at(0) * coefs.at(n).at(1));
                if (fabs(det) < eps)
                    continue;
                double x = -(coefs.at(m).at(1) * coefs.at(n).at(2) - coefs.at(m).at(2) * coefs.at(n).at(1)) / det;
                double y = (coefs.at(m).at(0) * coefs.at(n).at(2) - coefs.at(m).at(2) * coefs.at(n).at(0)) / det;
                bool is_local_voronoi_point = true;
                for ( unsigned int k = 0; k < coefs.size(); k++ )
                    if ( k != m && k != n )
                    {
                        if ( x * coefs.at(k).at(0) + y * coefs.at(k).at(1) + coefs.at(k).at(2) > eps )
                        {
                            is_local_voronoi_point = false;
                            break;
                        }
                    }
                if ( is_local_voronoi_point && vector2d(x, y) != get_point(i) )
                    if (std::find(voronoi_local_points.begin(), voronoi_local_points.end(), vector2d(x, y) - get_point(i)) == voronoi_local_points.end())
                        voronoi_local_points.push_back(vector2d(x, y) - get_point(i));
            }
        // sorting voronoi_local_points by angle with {-1, 0}
        std::sort(voronoi_local_points.begin(), voronoi_local_points.end(),
                  [] (const vector2d & p, const vector2d & q) -> bool
                  {
                        return p.angle() < q.angle();
                  });
        // calculating area and points and polygons of the voronoi diagram
        // note: need to process corner points (especially right upper: n = nx+ny) specially
        double area = 0;
        int voronn = voronp.size();
        voronp.push_back(get_point(i));
        voront.push_back({});
        for ( unsigned int j = 0; j < voronoi_local_points.size() - 1; j++ )
        {
            area += fabs( Vec(voronoi_local_points.at(j), voronoi_local_points.at(j+1)) )/2;
            voronp.push_back(get_point(i) + voronoi_local_points.at(j));
            voront.back().push_back(voronp.size() - 1);
        }
        voronp.push_back(voronoi_local_points.back() + get_point(i));
        voront.back().push_back(voronp.size() - 1);
        if ( !is_corner(i) )
        {
            area += fabs( Vec(voronoi_local_points.back(), voronoi_local_points.front()) )/2;
        }
        else if ( i == nx + ny )
        {
            area += fabs( Vec(voronoi_local_points.back(), voronoi_local_points.front()) )/2;
            area -= fabs( Vec(voronoi_local_points.back(), voronoi_local_points.at(voronoi_local_points.size()-2)) )/2;
            voront.back().insert(voront.back().end() - 1, voronn);
        }
        else
        {
            voront.back().push_back(voronn);
        }

        voronoi_areas.at(i) = area;
        coefs.clear();
        voronoi_local_points.clear();
    }

    /////////////////////////////////
    // Saving voronoi diagram
    vtk_file << "DATASET POLYDATA\nPOINTS " << voronp.size() << " float\n";
    for ( auto vp : voronp )
        vtk_file << vp.x << " " << vp.y << " "  << 0.0 << "\n";
    int vtsum = 0;
    for ( auto vt : voront ) vtsum += vt.size()+1;
    vtk_file << "\nPOLYGONS " << voront.size() << " " << vtsum << "\n";
    for ( auto vt : voront )
    {
        vtk_file << vt.size() << " ";
        for (auto el : vt) vtk_file << " " << el;
        vtk_file << "\n";
    }
    vtk_file << "\nPOINT_DATA " << voronp.size() << "\n";;
    vtk_file << "SCALARS voronz FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( unsigned int i = 0; i < voronp.size(); i++ )
    {
        vtk_file <<  0 << "\n";
    }
}


int mesh_2d::get_number_of_points()
{
	return points.size();
}

int mesh_2d::get_number_of_main_points()
{
    return points.size() / (N*N);
}

int mesh_2d::get_number_of_contour_points()
{
    return cont_points.size();
}

vector2d mesh_2d::get_point(int n)
{
    return points.at(n);
}

int mesh_2d::get_number_of_elements()
{
	return elements.size();
}

int mesh_2d::get_triangle_point_num(int n, int k)
{
    return elements.at(n).at(k);
}

vector2d mesh_2d::get_triangle_point(int n, int k)
{
	return get_point(get_triangle_point_num(n,k));
}

bool mesh_2d::is_inside(vector2d p)
{
    if (p.x < size_x0+eps || p.x > size_x1-eps || p.y < size_y0+eps || p.y > size_y1-eps)
		return false;
	else return true;
}

void mesh_2d::make_inside_continuous(vector2d & p)
{
    if (p.x < size_x0)
        p.x += size_x1-size_x0;
    else if (p.x > size_x1)
        p.x -= size_x1-size_x0;
    if (p.y < size_y0)
        p.y += size_y1 - size_y0;
    else if (p.y > size_y1)
        p.y -= size_y1 - size_y0;
}

void mesh_2d::make_inside_symmetric(vector2d & p)
{
    if (p.x < size_x0)
        p.x += 2*(size_x0-p.x);
    else if (p.x > size_x1)
        p.x += 2*(size_x1-p.x);
    if (p.y < size_y0)
        p.y += 2*(size_y0-p.y);
    else if (p.y > size_y1)
        p.y += 2*(size_y1-p.y);
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

//    vector2d p1 = get_triangle_point(n, 0);
//    vector2d p2 = get_triangle_point(n, 1);
//    vector2d p3 = get_triangle_point(n, 2);

//    double alpha = ((p2.y - p3.y)*(p.x - p3.x) + (p3.x - p2.x)*(p.y - p3.y)) /
//            ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y));
//    double beta = ((p3.y - p1.y)*(p.x - p3.x) + (p1.x - p3.x)*(p.y - p3.y)) /
//           ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y));
//    double gamma = 1.0f - alpha - beta;

//    double eps = 1e-4*h;

//    return alpha >= -eps && alpha <= 1.0+eps && beta >= -eps && beta <= 1.0+eps && gamma <= 1.0+eps;
}

vector2d mesh_2d::get_center_point(int element_num)
{
    return (points[elements[element_num][0]] +
            points[elements[element_num][1]] +
            points[elements[element_num][2]])/3.0;
}


bool mesh_2d::is_corner(int n)
{
	return (n == 0) || (n == nx) || (n == nx + ny) || (n == 2*nx + ny);
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
    {
        std::cout << "Error in 'get_opposite_point_num'!\n";
        std::exit(1);
    }
}

