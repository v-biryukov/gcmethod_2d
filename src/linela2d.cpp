///////////////////////////////////////////////////////////
//////
//////	file: linela2d.cpp
//////  class for solving linear elastisity equation in 2d using
//////  grid-characteristic method.
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////

#include "linela2d.h"
#include <sstream>

point_data linela2d::rotate(point_data & origin, int sign)
{
    double rxx = directions[2].x;
    double rxy = sign*directions[2].y;
    double ryx = sign*directions[4].x;
    double ryy = directions[4].y;
    point_data result;
    result.vx = rxx * origin.vx + rxy * origin.vy;
    result.vy = ryx * origin.vx + ryy * origin.vy;
    result.sxx = rxx*rxx*origin.sxx + 2*rxx*rxy*origin.sxy + rxy*rxy*origin.syy;
    result.sxy = rxx*ryx*origin.sxx + (rxx*ryy+rxy*ryx)*origin.sxy + rxy*ryy*origin.syy;
    result.syy = ryx*ryx*origin.sxx + 2*ryx*ryy*origin.sxy + ryy*ryy*origin.syy;
    return result;
}

riemann_data linela2d::get_riemann_inv_X(int pn)
{
    point_data pd;
    if (is_axes_random)
        pd = rotate(data[pn], 1);
    else
        pd = data[pn];
    riemann_data rd;

    int en1 = eldata_X[pn].el[0];
    int en2 = eldata_X[pn].el[1];
    //double c10 = (c1[en1] + c1[en2]) / 2.0;
    //double c30 = (c3(en1) + c3(en2)) / 2.0;

    rd.w[0] =  0;//-c30/c10 * pd.sxx   +     pd.syy;

    rd.w[1] = -0.5 * c2[en1]*rho[en1]  * pd.vy    +   0.5                     * pd.sxy;

    rd.w[2] = +0.5 * c2[en2]*rho[en2]  * pd.vy    +   0.5                     * pd.sxy;

    rd.w[3] = - 0.5 * rho[en1] * c3(en1) * pd.vx    +   0.5 * c3(en1) / c1[en1] * pd.sxx;

    rd.w[4] = + 0.5 * rho[en2] * c3(en2) * pd.vx    +   0.5 * c3(en2) / c1[en2] * pd.sxx;

    return rd;
}

void linela2d::set_point_data_X(int pn, riemann_data & rd)
{
    int en1 = eldata_X[pn].el[0];
    int en2 = eldata_X[pn].el[1];
    point_data pd;

    pd.vx = - 1.0/(rho[en1] * c3(en1))  * rd.w[3] + 1.0/(rho[en2] * c3(en2))   * rd.w[4];

    pd.vy = - 1.0 / (rho[en1] * c2[en1])   * rd.w[1] + 1.0 / (rho[en2] * c2[en2])  * rd.w[2];

    pd.sxx = +c1[en1]/c3(en1)   * rd.w[3] + c1[en2]/c3(en2)   * rd.w[4];

    pd.sxy =    rd.w[1] +   rd.w[2];

    pd.syy =    rd.w[0] +   rd.w[3] + rd.w[4];

    if (is_axes_random)
        pd = rotate(pd, -1);

    data[pn] += pd;
}

riemann_data linela2d::get_riemann_inv_Y(int pn)
{
    point_data pd;
    if (is_axes_random)
        pd = rotate(data[pn], 1);
    else
        pd = data[pn];
    riemann_data rd;

    int en1 = eldata_Y[pn].el[0];
    int en2 = eldata_Y[pn].el[1];
    //double c10 = (c1[en1] + c1[en2]) / 2.0;
    //double c30 = (c3(en1) + c3(en2)) / 2.0;

    rd.w[0] =  0.0;//1.0                    * pd.sxx   - c30 / c10 * pd.syy;

    rd.w[1] = -0.5 * c2[en1]*rho[en1] * pd.vx    +   0.5     * pd.sxy;

    rd.w[2] = +0.5 * c2[en2]*rho[en2] * pd.vx    +   0.5     * pd.sxy;

    rd.w[3] = -0.5 * c1[en1]*rho[en1] * pd.vy    +   0.5     * pd.syy;

    rd.w[4] = +0.5 * c1[en2]*rho[en2] * pd.vy    +   0.5     * pd.syy;
    return rd;
}

void linela2d::set_point_data_Y(int pn, riemann_data & rd)
{
    int en1 = eldata_Y[pn].el[0];
    int en2 = eldata_Y[pn].el[1];
    point_data pd;

    pd.vx = - 1.0 / (c2[en1] * rho[en1])  * rd.w[1] + 1.0 / (c2[en2] * rho[en2])  * rd.w[2];

    pd.vy = - 1.0 / (c1[en1] * rho[en1])  * rd.w[3] + 1.0 / (c1[en2] * rho[en2])  * rd.w[4];

    pd.sxx = rd.w[0] + c3(en1)/c1[en1]  * rd.w[3] + c3(en2)/c1[en2]  * rd.w[4];

    pd.sxy =                    rd.w[1] +                   rd.w[2];

    pd.syy =             rd.w[3] + rd.w[4];


    if (is_axes_random)
        pd = rotate(pd, -1);

    data[pn] += pd;
}

std::vector<double> linela2d::get_lambda_X(int point_n)
{
    std::vector<double> lambda(5);
    lambda[0] = 0.0;
    lambda[1] = - c2[eldata_X[point_n].el[0]];
    lambda[2] = + c2[eldata_X[point_n].el[1]];
    lambda[3] = - c1[eldata_X[point_n].el[0]];
    lambda[4] = + c1[eldata_X[point_n].el[1]];
    return lambda;
}

std::vector<double> linela2d::get_lambda_Y(int point_n)
{
    std::vector<double> lambda(5);
    lambda[0] = 0.0;
    lambda[1] = - c2[eldata_Y[point_n].el[0]];
    lambda[2] = + c2[eldata_Y[point_n].el[1]];
    lambda[3] = - c1[eldata_Y[point_n].el[0]];
    lambda[4] = + c1[eldata_Y[point_n].el[1]];
    return lambda;
}

void linela2d::calculate_point_elements()
{
    for (int i = 0; i < mesh->get_number_of_points(); i++)
        for (int k = 1; k < 5; k++)
        {
            vector2d p = mesh->get_point(i) + lambda_hint * tau * directions[k];
            if ( !mesh->is_inside(p) ) mesh->make_inside_vector(p);
            for (int j = 0; j < mesh->triangles[i].size(); j++)
            {
                int n = mesh->triangles[i][j];
                if (mesh->is_inside(p, n))
                {
                    if (k == 1)
                        eldata_X[i].el[0] = n;
                    else if (k==2)
                        eldata_X[i].el[1] = n;
                    else if (k==3)
                        eldata_Y[i].el[0] = n;
                    else if (k==4)
                        eldata_Y[i].el[1] = n;
                }
            }
        }
}

void linela2d::compute_rdata_X()
{
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        rdata[i] = get_riemann_inv_X(i);
    }
}

void linela2d::compute_rdata_Y()
{
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        rdata[i] = get_riemann_inv_Y(i);
    }
}

void linela2d::step_X()
{
    compute_rdata_X();
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        std::vector<double> lambda = get_lambda_X(i);
        riemann_data temp_rd;
        temp_rd.w[0] = 0;
        for (int k = 1; k < 5; k++)
        {
            vector2d p = mesh->get_point(i) + lambda[k] * tau * directions[2];
            if ( !mesh->is_inside(p) ) mesh->make_inside_vector(p);
            temp_rd.w[k] = approximate(p, eldata_X[i].el[(k-1)%2] , k) - rdata[i].w[k];
        }
        set_point_data_X(i, temp_rd);
    }
}

void linela2d::step_Y()
{
    compute_rdata_Y();
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        std::vector<double> lambda = get_lambda_Y(i);
        riemann_data temp_rd;
        temp_rd.w[0] = 0;
        for (int k = 1; k < 5; k++)
        {
            vector2d p = mesh->get_point(i) + lambda[k] * tau * directions[4];
            if ( !mesh->is_inside(p) ) mesh->make_inside_vector(p);
            temp_rd.w[k] = approximate(p, eldata_Y[i].el[(k-1)%2] , k) - rdata[i].w[k];
        }
        set_point_data_Y(i, temp_rd);
    }
}

void linela2d::step()
{
    if (is_axes_random)
    {
        set_directions(((double) rand() / (RAND_MAX)) * 2.0 * M_PI);
        calculate_point_elements();
    }
    step_X();
    step_Y();
}

void linela2d::calculate()
{

    for (int i = 0; i < number_of_steps; i++)
    {
        if (i % saving_frequency == 0)
        {
            std::stringstream ss;
            ss << i / saving_frequency;
            save_to_vtk("out/linela_" + ss.str() + ".vtk");
            std::cout << i << " steps from " << number_of_steps << std::endl;
        }
        step();
    }
}

double linela2d::approximate(vector2d r, int tn, int k)
{
    vector2d ra = mesh->get_triangle_point(tn, 0);
    vector2d rb = mesh->get_triangle_point(tn, 1);
    vector2d rc = mesh->get_triangle_point(tn, 2);
    double sa = Vec(rc - rb, r - rb)/2;
    double sb = Vec(ra - rc, r - rc)/2;
    double sc = Vec(rb - ra, r - ra)/2;
    double s = sa + sb + sc;
    sa /= s; sb /= s; sc /= s;
    double v = 0;
    if ( N == 1 )
    {
        v += sa * rdata[mesh->elements[tn][0]].w[k];
        v += sb * rdata[mesh->elements[tn][1]].w[k];
        v += sc * rdata[mesh->elements[tn][2]].w[k];

    }
    else if ( N == 2 )
    {
        v += sa*(2*sa-1) * rdata[mesh->elements[tn][0]].w[k];
        v += sb*(2*sb-1) * rdata[mesh->elements[tn][1]].w[k];
        v += sc*(2*sc-1) * rdata[mesh->elements[tn][2]].w[k];
        v += 4*sa*sb     * rdata[mesh->elements[tn][3]].w[k];
        v += 4*sb*sc     * rdata[mesh->elements[tn][4]].w[k];
        v += 4*sc*sa     * rdata[mesh->elements[tn][5]].w[k];
    }
    else if ( N == 3 )
    {
        v += sa*(3*sa-1)*(3*sa-2)/2 * rdata[mesh->elements[tn][0]].w[k];
        v += sb*(3*sb-1)*(3*sb-2)/2 * rdata[mesh->elements[tn][1]].w[k];
        v += sc*(3*sc-1)*(3*sc-2)/2 * rdata[mesh->elements[tn][2]].w[k];
        v += 9*sa*(3*sa-1)*sb/2     * rdata[mesh->elements[tn][3]].w[k];
        v += 9*sb*(3*sb-1)*sa/2     * rdata[mesh->elements[tn][4]].w[k];
        v += 9*sb*(3*sb-1)*sc/2     * rdata[mesh->elements[tn][5]].w[k];
        v += 9*sc*(3*sc-1)*sb/2     * rdata[mesh->elements[tn][6]].w[k];
        v += 9*sc*(3*sc-1)*sa/2     * rdata[mesh->elements[tn][7]].w[k];
        v += 9*sa*(3*sa-1)*sc/2     * rdata[mesh->elements[tn][8]].w[k];
        v += 27*sa*sb*sc            * rdata[mesh->elements[tn][9]].w[k];
    }
    else if ( N == 4 )
    {
        v += sa*(4*sa-1)*(2*sa-1)*(4*sa-3)/3 * rdata[mesh->elements[tn][0]].w[k];
        v += sb*(4*sb-1)*(2*sb-1)*(4*sb-3)/3 * rdata[mesh->elements[tn][1]].w[k];
        v += sc*(4*sc-1)*(2*sc-1)*(4*sc-3)/3 * rdata[mesh->elements[tn][2]].w[k];
        v += 16*sa*(4*sa-1)*(2*sa-1)*sb/3    * rdata[mesh->elements[tn][3]].w[k];
        v += 4*sa*(4*sa-1)*sb*(4*sb-1)       * rdata[mesh->elements[tn][4]].w[k];
        v += 16*sb*(4*sb-1)*(2*sb-1)*sa/3    * rdata[mesh->elements[tn][5]].w[k];
        v += 16*sb*(4*sb-1)*(2*sb-1)*sc/3    * rdata[mesh->elements[tn][6]].w[k];
        v += 4*sb*(4*sb-1)*sc*(4*sc-1)       * rdata[mesh->elements[tn][7]].w[k];
        v += 16*sc*(4*sc-1)*(2*sc-1)*sb/3    * rdata[mesh->elements[tn][8]].w[k];
        v += 16*sc*(4*sc-1)*(2*sc-1)*sa/3    * rdata[mesh->elements[tn][9]].w[k];
        v += 4*sa*(4*sa-1)*sc*(4*sc-1)       * rdata[mesh->elements[tn][10]].w[k];
        v += 16*sa*(4*sa-1)*(2*sa-1)*sc/3    * rdata[mesh->elements[tn][11]].w[k];
        v += 32*sa*(4*sa-1)*sb*sc            * rdata[mesh->elements[tn][12]].w[k];
        v += 32*sb*(4*sb-1)*sc*sa            * rdata[mesh->elements[tn][13]].w[k];
        v += 32*sc*(4*sc-1)*sa*sb            * rdata[mesh->elements[tn][14]].w[k];
    }
    return v;
}

void linela2d::set_directions(double angle)
{
    directions[2] = vector2d(cos(angle), sin(angle));
    directions[1] = - directions[2];
    directions[4] = vector2d(-sin(angle), cos(angle));
    directions[3] = - directions[4];
}

void linela2d::init()
{
    c1.resize(mesh->get_number_of_triangles());
    c2.resize(mesh->get_number_of_triangles());
    rho.resize(mesh->get_number_of_triangles());
    data.resize(mesh->get_number_of_points());
    //data_prev.resize(mesh->get_number_of_points());
    rdata.resize(mesh->get_number_of_points());
    eldata_X.resize(mesh->get_number_of_points());
    eldata_Y.resize(mesh->get_number_of_points());
    lambda_hint = mesh->get_max_altitude() / N / 2.0;

    directions.reserve(5);
    directions.push_back(vector2d(0, 0));
    directions.push_back(vector2d(-1, 0));
    directions.push_back(vector2d(1, 0));
    directions.push_back(vector2d(0, -1));
    directions.push_back(vector2d(0, 1));

    calculate_point_elements();
}

void linela2d::read_from_file(std::string path)
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

    tau = pt.get<double>("Method.tau");
    number_of_steps = pt.get<int>("Method.number_of_steps");
    N = pt.get<int>("Method.N");

    string is_monotonic_str = pt.get<string>("Method.is_monotonic");
    is_monotonic = (is_monotonic_str == "true" || is_monotonic_str == "TRUE" || is_monotonic_str == "True");

    string is_axes_random_str = pt.get<string>("Method.is_axes_random");
    is_axes_random = (is_axes_random_str == "true" || is_axes_random_str == "TRUE" || is_axes_random_str == "True");

    //string use_precalculations_str = pt.get<string>("Method.use_precalculations");
    //use_precalculations = (use_precalculations_str == "true" || use_precalculations_str == "TRUE" || use_precalculations_str == "True");

    //string save_only_main_points_str = pt.get<string>("Method.save_only_main_points");
    //save_only_main_points = (save_only_main_points_str == "true" || save_only_main_points_str == "TRUE" || save_only_main_points_str == "True");

    saving_frequency = pt.get<int>("Method.saving_frequency");
}

void linela2d::save_to_vtk(std::string name)
{
    std::ofstream vtk_file(name.c_str());
    vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
    vtk_file << "DATASET POLYDATA\nPOINTS " << mesh->get_number_of_points() << " float\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
        vtk_file << mesh->get_point(i).x << " " << mesh->get_point(i).y << " "  << 0.0 << "\n";
    vtk_file << "\nPOLYGONS " << mesh->get_number_of_triangles() << " " << mesh->get_number_of_triangles()*4 << "\n";
    for (int i = 0; i < mesh->get_number_of_triangles(); i++)
        vtk_file << 3 << " " << mesh->get_triangle_point_num(i,0) << " " << mesh->get_triangle_point_num(i,1) << " " << mesh->get_triangle_point_num(i,2) << "\n";
    vtk_file << "\nPOINT_DATA " << mesh->get_number_of_points() << "\n";;
    vtk_file << "VECTORS v FLOAT\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        vtk_file <<  data[i].vx << " " << data[i].vy << " " << 0.0 << "\n";
    }
    vtk_file << "SCALARS t11 FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        vtk_file <<  data[i].sxx << "\n";
    }
    vtk_file << "SCALARS t12 FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        vtk_file <<  data[i].sxy << "\n";
    }
    vtk_file << "SCALARS t22 FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        vtk_file <<  data[i].syy << "\n";
    }
}
