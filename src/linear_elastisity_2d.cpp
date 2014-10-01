///////////////////////////////////////////////////////////
//////
//////	file: linear_elastisity_2d.cpp
//////  class for solving linear elastisity equation in 2d using
//////  grid-characteristic method.
//////  author: Biryukov Vladimir, biryukov.vova@gmail.com
//////  MIPT, 2014
//////
///////////////////////////////////////////////////////////

#include "linear_elastisity_2d.h"



void DrawProgressBar(int len, double percent) {
  std::cout << "\x1B[2K"; // Erase the entire current line.
  std::cout << "\x1B[0E"; // Move to the beginning of the current line.
  std::string progress;
  for (int i = 0; i < len; ++i) {
    if (i < static_cast<int>(len * percent)) {
      progress += "=";
    } else {
      progress += " ";
    }
  }
  std::cout << "[" << progress << "] " << (static_cast<int>(100 * percent)) << "%";
  flush(std::cout); // Required.
}


void linear_elastisity_2d::initial_conditions(vector2d & v, tensor2d & T, double x, double y)
{
    set_S_wave(v, T, x, y, vector2d(0, 0), vector2d(0, -1), 2, 1);
    /*
    if ( x > -1.5 && x < 1.5)
        v.x = 2;
    else
        v.x = 0;
    v.y = 0;

    if ( x*x + y*y < 1 )
    {
        v.x = 2;
        v.y = 0;
    }
    else
    {
        v.x = 0;
        v.y = 0;
    }

    T.t00 = 0;
    T.t01 = 0;
    T.t10 = T.t01;
    T.t11 = 0;*/
}

void linear_elastisity_2d::initial_conditions(int j, std::vector<double> & w, double x, double y)
{
    vector2d v;
    tensor2d T;
    initial_conditions(v, T, x, y);
    get_w(j, w, v, T);
}
double linear_elastisity_2d::initial_conditions(int j, int k, double x, double y)
{
    vector2d v;
    tensor2d T;
    initial_conditions(v, T, x, y);
    return get_w(j, k, v, T);
}

void linear_elastisity_2d::linear_elastisity_2d::init()
{
    ns.resize(2, std::vector<vector2d>(2));
    ls.resize(2);
    ns[0][0] = xis[0]/Magnitude(xis[0]);
    ns[1][0] = xis[1]/Magnitude(xis[1]);
    for ( int i = 0; i < 2; i++ )
    {
        ls[i] = Magnitude(xis[i]);
        ns[i][1] = vector2d(-ns[i][0].y, ns[i][0].x);
    }

    N0.resize(2, std::vector<tensor2d>(2));
    for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ )
            N0[i][j] = ((ns[i][0]^ns[j][0]) + (ns[j][0]^ns[i][0]))/2;
    N1.resize(2, std::vector<tensor2d>(2));
    for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ )
            N1[i][j] = ((ns[i][0]^ns[j][0]) + (ns[j][0]^ns[i][0]))/2;

    Lambda.resize(5);
    Lambda.at(0) =  c1 ;
    Lambda.at(1) = -c1 ;
    Lambda.at(2) =  c2 ;
    Lambda.at(3) = -c2 ;
    Lambda.at(4) =   0 ;

    values_v0.resize(mesh->get_number_of_points());
    values_T0.resize(mesh->get_number_of_points());
    values_v1.resize(mesh->get_number_of_points());
    values_T1.resize(mesh->get_number_of_points());
    set_initial_conditions();
    cur_step = 0;
}

double linear_elastisity_2d::c3()
{
    return c1*(1 - 2 * c2 * c2 / c1 / c1);
}

void linear_elastisity_2d::set_P_wave(vector2d & v, tensor2d & T, double x, double y, vector2d c, vector2d dir, double w, double mag)
{
    dir /= Magnitude(dir);
    vector2d r = vector2d(x, y);
    tensor2d I = tensor2d(1, 0, 0, 1);
    double d = (r-c)*dir;
    if ( abs(d) < w/2 )
    {
        mag *= cos(3.14159265*d/w);
        v = dir * mag;
        T = - mag * rho * c3() * I;
        T -= mag * rho * (c1 - c3()) * tensor2d(dir.x*dir.x, dir.x*dir.y, dir.x*dir.y, dir.y*dir.y);
    }
    else
    {
        v = vector2d(0, 0);
        T = tensor2d(0, 0, 0, 0);
    }
}

void linear_elastisity_2d::set_S_wave(vector2d & v, tensor2d & T, double x, double y, vector2d c, vector2d dir, double w, double mag)
{
    vector2d r = vector2d(x, y);
    vector2d ort_dir = vector2d(-dir.y, dir.x);
    double d = (r-c)*dir;
    if ( abs(d) < w/2 )
    {
        mag *= cos(3.14159265*d/w);
        v = mag * ort_dir;
        T = - mag * rho * c2 * tensor2d(dir.x*v.x, dir.x*v.y, dir.y*v.x, dir.y*v.y);
    }
    else
    {
        v = vector2d(0, 0);
        T = tensor2d(0, 0, 0, 0);
    }
}

void linear_elastisity_2d::set_initial_conditions()
{
    vector2d v;
    tensor2d T;

    for ( unsigned int j = 0; j < mesh->points.size(); j++ )
    {
        initial_conditions(v, T, mesh->points.at(j).x, mesh->points.at(j).y);
        values_v0.at(j) = v;
        values_T0.at(j) = T;
    }
}

void linear_elastisity_2d::step()
{
    vector2d p;
    std::vector<std::vector<double> >  w;
    w.resize(2, std::vector<double>(5));
    std::vector<double> w_prev((N+1)*(N+2)/2);
    for ( int j = 0; j < 2; j++ )
    {
        for ( int i = 0; i < mesh->get_number_of_points(); i++ )
        {
            for ( int k = 0; k < 4; k++ )
            {
                int tn = 0;
                p = mesh->get_point(i) - Lambda.at(k)*xis.at(j)/Magnitude(xis.at(j))*tau;
                if ( !mesh->is_inside(p) ) mesh->make_inside_vector(p);
                for ( int n : mesh->triangles.at(i) )
                    if (mesh->is_inside(p, n))
                    {
                        tn = n;
                        for ( unsigned int m = 0; m < mesh->elements.at(n).size(); m++ )
                        {
                            int pn = mesh->elements.at(n).at(m);
                            w_prev.at(m) = get_w(j, k, values_v0.at(pn), values_T0.at(pn));
                        }
                        break;
                    }
                w.at(j).at(k) = approximate(w_prev, p, tn);
            }
            w.at(j).at(4) = initial_conditions(j, 4, mesh->points.at(i).x, mesh->points.at(i).y);
            get_vT(j, w.at(j), values_v1.at(i), values_T1.at(i));
        }
        values_v0.swap(values_v1);
        values_T0.swap(values_T1);
    }
    cur_step++;
}

double linear_elastisity_2d::approximate(std::vector<double> & w_prev, vector2d r, int tn)
{
	vector2d ra = mesh->get_triangle_point(tn, 0);
	vector2d rb = mesh->get_triangle_point(tn, 1);
	vector2d rc = mesh->get_triangle_point(tn, 2);
	double sa = Vec(rc - rb, r - rb)/2;
	double sb = Vec(ra - rc, r - rc)/2;
	double sc = Vec(rb - ra, r - ra)/2;
	double s = sa + sb + sc;
	sa /= s; sb /= s; sc /= s;
    std::vector<double> w;
	double v = 0;
	if ( N == 1 )
	{
        v += sa * w_prev.at(0);
        v += sb * w_prev.at(1);
        v += sc * w_prev.at(2);
    }
    else if ( N == 2 )
    {
        v += sa*(2*sa-1) * w_prev.at(0);
        v += sb*(2*sb-1) * w_prev.at(1);
        v += sc*(2*sc-1) * w_prev.at(2);
        v += 4*sa*sb * w_prev.at(3);
        v += 4*sb*sc * w_prev.at(4);
        v += 4*sc*sa * w_prev.at(5);
    }
    else if ( N == 3 )
    {
        v += sa*(3*sa-1)*(3*sa-2)/2 * w_prev.at(0);
        v += sb*(3*sb-1)*(3*sb-2)/2 * w_prev.at(1);
        v += sc*(3*sc-1)*(3*sc-2)/2 * w_prev.at(2);
        v += 9*sa*(3*sa-1)*sb/2     * w_prev.at(3);
        v += 9*sb*(3*sb-1)*sa/2     * w_prev.at(4);
        v += 9*sb*(3*sb-1)*sc/2     * w_prev.at(5);
        v += 9*sc*(3*sc-1)*sb/2     * w_prev.at(6);
        v += 9*sc*(3*sc-1)*sa/2     * w_prev.at(7);
        v += 9*sa*(3*sa-1)*sc/2     * w_prev.at(8);
        v += 27*sa*sb*sc            * w_prev.at(9);
    }
    else if ( N == 4 )
    {
        v += sa*(4*sa-1)*(2*sa-1)*(4*sa-3)/3 * w_prev.at(0);
        v += sb*(4*sb-1)*(2*sb-1)*(4*sb-3)/3 * w_prev.at(1);
        v += sc*(4*sc-1)*(2*sc-1)*(4*sc-3)/3 * w_prev.at(2);
        v += 16*sa*(4*sa-1)*(2*sa-1)*sb/3 * w_prev.at(3);
        v += 4*sa*(4*sa-1)*sb*(4*sb-1)    * w_prev.at(4);
        v += 16*sb*(4*sb-1)*(2*sb-1)*sa/3 * w_prev.at(5);
        v += 16*sb*(4*sb-1)*(2*sb-1)*sc/3 * w_prev.at(6);
        v += 4*sb*(4*sb-1)*sc*(4*sc-1)    * w_prev.at(7);
        v += 16*sc*(4*sc-1)*(2*sc-1)*sb/3 * w_prev.at(8);
        v += 16*sc*(4*sc-1)*(2*sc-1)*sa/3 * w_prev.at(9);
        v += 4*sa*(4*sa-1)*sc*(4*sc-1)    * w_prev.at(10);
        v += 16*sa*(4*sa-1)*(2*sa-1)*sc/3 * w_prev.at(11);
        v += 32*sa*(4*sa-1)*sb*sc         * w_prev.at(12);
        v += 32*sb*(4*sb-1)*sc*sa         * w_prev.at(13);
        v += 32*sc*(4*sc-1)*sa*sb         * w_prev.at(14);
    }
    if (is_monotonic && min_max_check(v, w_prev))
        return sa * w_prev.at(0) + sb * w_prev.at(1) + sc * w_prev.at(2);
    else
        return v;
}

bool linear_elastisity_2d::min_max_check(double z, std::vector<double> w_prev)
{
    return   (  z > w_prev.at(0) + eps_z
             && z > w_prev.at(1) + eps_z
             && z > w_prev.at(2) + eps_z )
           ||(  z < w_prev.at(0) - eps_z
             && z < w_prev.at(1) - eps_z
             && z < w_prev.at(2) - eps_z );
}


void linear_elastisity_2d::calculate()
{
    for ( int i = 0; i < number_of_steps; i++ )
    {
        DrawProgressBar(40, i*1.0/number_of_steps);
        if ( i % saving_frequency == 0 )
            save_to_vtk("out/linela_" + std::to_string(i/saving_frequency) + ".vtk");
        step();
    }
}

void linear_elastisity_2d::save_to_vtk(std::string name)
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
        vtk_file <<  values_v0.at(i).x << " " << values_v0.at(i).y << " " << 0.0 << "\n";
	}
	vtk_file << "SCALARS t11 FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		vtk_file <<  values_T0.at(i).t00 << "\n";
    }
    vtk_file << "SCALARS t12 FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		vtk_file <<  values_T0.at(i).t01 << "\n";
    }
    vtk_file << "SCALARS t22 FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		vtk_file <<  values_T0.at(i).t11 << "\n";
    }
}

void linear_elastisity_2d::get_w(int j, std::vector<double> & w, vector2d v, tensor2d T)
{
    w.resize(5);
    if (j == 0)
    {
        w[0] = ns[0][0] * v - 1.0/c1/rho * N0[0][0] * T;
        w[1] = ns[0][0] * v + 1.0/c1/rho * N0[0][0] * T;
        w[2] = ns[0][1] * v - 1.0/c2/rho * N0[0][1] * T;
        w[3] = ns[0][1] * v + 1.0/c2/rho * N0[0][1] * T;
        w[4] = (N0[1][1] - (1-2*c2*c2/c1/c1)*N0[0][0]) * T;
    }
    else if (j == 1)
    {
        w[0] = ns[1][0] * v - 1.0/c1/rho * N1[0][0] * T;
        w[1] = ns[1][0] * v + 1.0/c1/rho * N1[0][0] * T;
        w[2] = ns[1][1] * v - 1.0/c2/rho * N1[0][1] * T;
        w[3] = ns[1][1] * v + 1.0/c2/rho * N1[0][1] * T;
        w[4] = (N1[1][1] - (1-2*c2*c2/c1/c1)*N1[0][0]) * T;
    }
}

double linear_elastisity_2d::get_w(int j, int k, vector2d v, tensor2d T)
{
    double w = 0;
    if (j == 0)
    {
        switch (k)
        {
            case 0 : w = ns[0][0] * v - 1.0/c1/rho * N0[0][0] * T; break;
            case 1 : w = ns[0][0] * v + 1.0/c1/rho * N0[0][0] * T; break;
            case 2 : w = ns[0][1] * v - 1.0/c2/rho * N0[0][1] * T; break;
            case 3 : w = ns[0][1] * v + 1.0/c2/rho * N0[0][1] * T; break;
            case 4 : w = (N0[1][1] - (1-2*c2*c2/c1/c1)*N0[0][0]) * T; break;
        }
    }
    else if (j == 1)
    {
        switch (k)
        {
            case 0 : w = ns[1][0] * v - 1.0/c1/rho * N1[0][0] * T; break;
            case 1 : w = ns[1][0] * v + 1.0/c1/rho * N1[0][0] * T; break;
            case 2 : w = ns[1][1] * v - 1.0/c2/rho * N1[0][1] * T; break;
            case 3 : w = ns[1][1] * v + 1.0/c2/rho * N1[0][1] * T; break;
            case 4 : w = (N1[1][1] - (1-2*c2*c2/c1/c1)*N1[0][0]) * T; break;
        }
    }
    return w;
}

void linear_elastisity_2d::get_vT (int j, std::vector<double> w, vector2d & v, tensor2d & T)
{
    tensor2d I = tensor2d(1, 0, 0, 1);
    if (j == 0)
    {
        v = ((w[0]+w[1])*ns[0][0] + (w[2]+w[3])*ns[0][1])/2;
        T = (rho*(w[1]-w[0])*((c1-c3())*N0[0][0] + c2*I) + 2*c2*rho*(w[3]-w[2])*N0[0][1] + 2*w[4]*(I-N0[0][0]))/2;
    }
    else if (j == 1)
    {
        v = ((w[0]+w[1])*ns[1][0] + (w[2]+w[3])*ns[1][1])/2;
        T = (rho*(w[1]-w[0])*((c1-c3())*N1[0][0] + c2*I) + 2*c2*rho*(w[3]-w[2])*N1[0][1] + 2*w[4]*(I-N1[0][0]))/2;
    }
}

void linear_elastisity_2d::read_from_file(std::string path)
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
    xis.resize(2);
    std::stringstream ss;
    ss << pt.get<std::string>("Method.axis1");
    ss >> xis[0].x >> xis[0].y;
    ss.str("");
    ss.clear();
    ss << pt.get<std::string>("Method.axis2");
    ss >> xis[1].x >> xis[1].y;
    tau = pt.get<double>("Method.tau");
    c1 = pt.get<double>("Equation.c1");
    c2 = pt.get<double>("Equation.c2");
    rho = pt.get<double>("Equation.rho");
    number_of_steps = pt.get<int>("Method.number_of_steps");
    N = pt.get<int>("Method.N");
    std::string is_monotonic_str = pt.get<std::string>("Method.is_monotonic");
    is_monotonic = (is_monotonic_str == "true" || is_monotonic_str == "TRUE" || is_monotonic_str == "True");
    std::string use_precalculations_str = pt.get<std::string>("Method.use_precalculations");
    use_precalculations = (use_precalculations_str == "true" || use_precalculations_str == "TRUE" || use_precalculations_str == "True");
    std::string save_only_main_points_str = pt.get<std::string>("Method.save_only_main_points");
    save_only_main_points = (save_only_main_points_str == "true" || save_only_main_points_str == "TRUE" || save_only_main_points_str == "True");
    saving_frequency = pt.get<int>("Method.saving_frequency");
}



