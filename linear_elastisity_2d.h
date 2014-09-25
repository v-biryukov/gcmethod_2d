#include "vector2.h"
#include "mesh_2d.h"
#include "gcmethod_2d.h"
#include "tensor2d.h"


class linear_elastisity_2d
{
	mesh_2d * mesh;
    std::vector<gcmethod_2d> equations;

    double c1, c2, rho;

    std::vector<vector2d> xis;
    std::vector<std::vector<vector2d> > ns;
    std::vector<double> ls;
    std::vector<double> Lambda;
    std::vector<std::vector<tensor2d> > N0;
    std::vector<std::vector<tensor2d> > N1;

    std::vector<double> w0;
    std::vector<double> w1;


    int N;
	bool is_monotonic;
    bool use_precalculations;
    bool save_only_main_points;
    int saving_frequency;
    double tau, h;
	int number_of_steps;
	vector2d direction;

    void init();
    void get_w(int j, std::vector<double> & w, vector2d v, tensor2d T);
    void get_vT (int j, std::vector<double> w, vector2d & v, tensor2d & T);
    void set_initial_conditions();
    void read_from_file(std::string path);
    void save_to_vtk(std::string name);
    void initial_conditions(vector2d & v, tensor2d & T, double x, double y)
    {
        v.x = 2*exp(-x*x/5 - y*y/6);
        v.y = 2*exp(-x*x/5 - y*y/6);
        T.t00 = 2*exp(-x*x/5 - y*y/6);
        T.t01 = 2*exp(-x*x/5 - y*y/6);
        T.t10 = T.t01;
        T.t11 = 2*exp(-x*x/5 - y*y/6);
    }

    void initial_conditions(int j, std::vector<double> & w, double x, double y)
    {
        vector2d v;
        tensor2d T;
        initial_conditions(v, T, x, y);
        get_w(j, w, v, T);
    }

public:
    void calculate();
    linear_elastisity_2d (std::string path, mesh_2d * m)
    {
        mesh = m;
        read_from_file(path);
        init();
    }

};

void linear_elastisity_2d::init()
{
    ns.resize(2, std::vector<vector2d>(2));
    ls.resize(2);
    for ( int i = 0; i < 2; i++ )
    {
        ns[i][0] = vector2d(-xis[i].y, xis[i].x);
        ls[i] = Magnitude(ns[i][0]);
        ns[i][0] /= ls[i];
        ns[i][1] = vector2d(-ns[i][1].y, ns[i][1].x);
    }

    N0.resize(2, std::vector<tensor2d>(2));
    for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ )
            N0[i][j] = ((ns[0][i]^ns[0][j]) + (ns[0][j]^ns[0][i]))/2;
    N1.resize(2, std::vector<tensor2d>(2));
    for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ )
            N1[i][j] = ((ns[1][i]^ns[1][j]) + (ns[1][j]^ns[1][i]))/2;

    Lambda.resize(5);
    Lambda.at(0) = c1;
    Lambda.at(1) = -c1;
    Lambda.at(2) = c2;
    Lambda.at(3) = -c2;
    Lambda.at(4) = 0;

    struct convection_solver_gc_settings c;
    c.tau = tau;
    c.number_of_steps = saving_frequency;
    c.N = N;
	c.is_monotonic = is_monotonic;
    c.use_precalculations = use_precalculations;
    c.save_only_main_points = true;
    c.saving_frequency = 0;

    for ( int i = 0; i < 8; i++ )
    {
        c.lambda = Lambda.at(i/2) * ls.at(i%2);
        c.direction = xis[i%2]/Magnitude(xis[i%2]);
        gcmethod_2d g = gcmethod_2d(c, mesh);
        equations.push_back(g);
    }
    set_initial_conditions();
}

void linear_elastisity_2d::set_initial_conditions()
{
    std::vector<double> w(5);

    for ( int j = 0; j < mesh->points.size(); j++ )
    {
        for ( int i = 0; i < 8; i++ )
        {
            initial_conditions(i%2, w, mesh->points.at(i).x, mesh->points.at(i).y);
            equations[i].values0[j] = w[i/2];
        }
    }
}

void linear_elastisity_2d::calculate()
{
    for ( int i = 0; i < number_of_steps/saving_frequency; i++ )
    {
        for ( int i = 0; i < 8; i++ )
        {
            equations[i].calculate();
        }
        save_to_vtk("out/linela_" + std::to_string(i) + ".vtk");
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
	vtk_file << "SCALARS vx FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
        initial_conditions(0, w0, mesh->points.at(i).x, mesh->points.at(i).y);
        initial_conditions(0, w1, mesh->points.at(i).x, mesh->points.at(i).y);
        for ( int i = 0; i < 4; i++ )
            w0[i] =  equations[2*i].values0.at(i);
        for ( int i = 0; i < 4; i++ )
            w1[i] = equations[2*i+1].values0.at(i);
        vector2d v0, v1;
        tensor2d T0, T1;
        get_vT(0, w0, v0, T0);
        get_vT(1, w1, v1, T1);
        vtk_file <<  (v0.x+v1.x)/2 << "\n";
	}
	vtk_file << "SCALARS vy FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
        initial_conditions(0, w0, mesh->points.at(i).x, mesh->points.at(i).y);
        initial_conditions(0, w1, mesh->points.at(i).x, mesh->points.at(i).y);
        for ( int i = 0; i < 4; i++ )
            w0[i] =  equations[2*i].values0.at(i);
        for ( int i = 0; i < 4; i++ )
            w1[i] = equations[2*i+1].values0.at(i);
        vector2d v0, v1;
        tensor2d T0, T1;
        get_vT(0, w0, v0, T0);
        get_vT(1, w1, v1, T1);
        vtk_file <<  (v0.y+v1.y)/2 << "\n";
    }
}

void linear_elastisity_2d::get_w(int j, std::vector<double> & w, vector2d v, tensor2d T)
{
    w.resize(5);
    if (j == 0)
    {
        w[0] = ns[0][0] * v - 1.0/c1/rho * N0[0][0] * T;
        w[1] = ns[0][0] * v + 1.0/c1/rho * N0[0][0] * T;
        w[2] = ns[0][1] * v - 1.0/c1/rho * N0[0][1] * T;
        w[3] = ns[0][1] * v + 1.0/c1/rho * N0[0][1] * T;
        w[4] = (N0[1][1] - (1-2*c2*c2/c1/c1)*N0[0][0]) * T;
    }
    else if (j == 1)
    {
        w[0] = ns[1][0] * v - 1.0/c1/rho * N1[0][0] * T;
        w[1] = ns[1][0] * v + 1.0/c1/rho * N1[0][0] * T;
        w[2] = ns[1][1] * v - 1.0/c1/rho * N1[0][1] * T;
        w[3] = ns[1][1] * v + 1.0/c1/rho * N1[0][1] * T;
        w[4] = (N1[1][1] - (1-2*c2*c2/c1/c1)*N1[0][0]) * T;
    }
}

void linear_elastisity_2d::get_vT (int j, std::vector<double> w, vector2d & v, tensor2d & T)
{
    tensor2d I = tensor2d(1, 0, 0, 1);
    if (j == 0)
    {
        v = (w[0]+w[1])*ns[0][0] + (w[2]+w[3])*ns[0][1];
        T = rho*(w[1]-w[0])*((c1-c2)*N0[0][0] + c2*I) + 2*c2*rho*(w[3]-w[2])*N0[0][1] + 2*w[4]*(I-N0[0][0]);
    }
    else if (j == 1)
    {
        v = (w[0]+w[1])*ns[1][0] + (w[2]+w[3])*ns[1][1];
        T = rho*(w[1]-w[0])*((c1-c2)*N1[0][0] + c2*I) + 2*c2*rho*(w[3]-w[2])*N1[0][1] + 2*w[4]*(I-N1[0][0]);
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
    std::string direction_str = pt.get<std::string>("Equation.direction");
    std::stringstream ss;
    ss << direction_str;
    ss >> direction.x >> direction.y;
    xis.resize(2);
    ss.str("");
    ss.clear();
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


