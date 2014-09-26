#include "vector2.h"
#include "mesh_2d.h"
#include "gcmethod_2d.h"
#include "tensor2d.h"


class linear_elastisity_2d
{
	mesh_2d * mesh;

	std::vector<vector2d> values_v0;
	std::vector<tensor2d> values_T0;
	std::vector<vector2d> values_v1;
	std::vector<tensor2d> values_T1;

    double c1, c2, rho;

    std::vector<vector2d> xis;
    std::vector<std::vector<vector2d> > ns;
    std::vector<double> ls;
    std::vector<double> Lambda;
    std::vector<std::vector<tensor2d> > N0;
    std::vector<std::vector<tensor2d> > N1;


    int N;
	bool is_monotonic;
    bool use_precalculations;
    bool save_only_main_points;
    int saving_frequency;
    double tau, h;
	int number_of_steps;

    void init();
    void get_w(int j, std::vector<double> & w, vector2d v, tensor2d T);
    void get_vT (int j, std::vector<double> w, vector2d & v, tensor2d & T);
    void set_initial_conditions();
    void read_from_file(std::string path);
    void save_to_vtk(std::string name);
    void initial_conditions(vector2d & v, tensor2d & T, double x, double y)
    {
        if ( x*x + y*y < 1 )
        //if ( x > -1 && x < 1)
            v.x = 2;
        else
            v.x = 0;
        v.y = 0;
        T.t00 = 0;
        T.t01 = 0;
        T.t10 = T.t01;
        T.t11 = 0;
    }

    void initial_conditions(int j, std::vector<double> & w, double x, double y)
    {
        vector2d v;
        tensor2d T;
        initial_conditions(v, T, x, y);
        get_w(j, w, v, T);
    }
    double initial_conditions(int j, int k, double x, double y)
    {
        vector2d v;
        tensor2d T;
        initial_conditions(v, T, x, y);
        return get_w(j, k, v, T);
    }

    void step();
    double approximate(std::vector<double> & w_prev, vector2d r, int tn);
    double get_w(int j, int k, vector2d v, tensor2d T);
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
    double xi_det = xis[0].x*xis[2].y - xis[0].y*xis[2].x;
    ns[0][0] = vector2d(xis[2].y, -xis[2].x) / xi_det;
    ns[1][0] = vector2d(-xis[0].y, xis[0].x) / xi_det;
    for ( int i = 0; i < 2; i++ )
    {
        ls[i] = Magnitude(ns[i][0]);
        ns[i][0] /= ls[i];
        ns[i][1] = vector2d(-ns[i][0].y, ns[i][0].x);
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
}

void linear_elastisity_2d::set_initial_conditions()
{
    vector2d v;
    tensor2d T;

    for ( int j = 0; j < mesh->points.size(); j++ )
    {
        initial_conditions(v, T, mesh->points.at(j).x, mesh->points.at(j).y);
        values_v0.at(j) = v;
        values_T0.at(j) = T;
    }
}

void linear_elastisity_2d::step()
{
    vector2d p;
    for ( int j = 0; j < 2; j++ )
    {
        for ( int i = 0; i < mesh->get_number_of_points(); i++ )
        {
            std::vector<std::vector<double> >  w;
            w.resize(2, std::vector<double>(5));
            for ( int k = 0; k < 4; k++ )
            {
                int tn;
                std::vector<double> w_prev((N+1)*(N+2)/2);
                p = mesh->get_point(i) - Lambda.at(k)*xis.at(k)*tau;
                if (p.x > -0.5 && p.x < 0.5)
                    int aaaa = 1;
                if ( !mesh->is_inside(p) ) mesh->make_inside_vector(p);
                for ( int n : mesh->triangles.at(i) )
                    if (mesh->is_inside(p, n))
                    {
                        tn = n;
                        for ( int m = 0; m < mesh->elements.at(n).size(); m++ )
                        {
                            int pn = mesh->elements.at(n).at(m);
                            w_prev.at(m) = get_w(j, k, values_v0.at(pn), values_T0.at(pn));
                        }
                        break;
                    }
                w.at(j).at(k) = approximate(w_prev, p, tn);
            }
            w.at(j).at(4) = initial_conditions(j, 4, mesh->points.at(i).x, mesh->points.at(i).y);
            vector2d v;
            tensor2d T;
            get_vT(j, w.at(j), v, T);
            values_v1.at(i) = v;
            values_T1.at(i) = T;
        }
        values_v0.swap(values_v1);
        values_T0.swap(values_T1);
    }
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
	//std::cout << w_prev.at(0) << " " << w_prev.at(1) << " " << w_prev.at(2) << "\n";
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
	return v;
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
	vtk_file << "SCALARS vx FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
        vtk_file <<  values_v0.at(i).x << "\n";
	}
	vtk_file << "SCALARS vy FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh->get_number_of_points(); i++ )
	{
		vtk_file <<  values_v0.at(i).y << "\n";
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
    double w;
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
        T = (rho*(w[1]-w[0])*((c1-c2)*N0[0][0] + c2*I) + 2*c2*rho*(w[3]-w[2])*N0[0][1] + 2*w[4]*(I-N0[0][0]))/2;
    }
    else if (j == 1)
    {
        v = ((w[0]+w[1])*ns[1][0] + (w[2]+w[3])*ns[1][1])/2;
        T = (rho*(w[1]-w[0])*((c1-c2)*N1[0][0] + c2*I) + 2*c2*rho*(w[3]-w[2])*N1[0][1] + 2*w[4]*(I-N1[0][0]))/2;
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
    xis.resize(4);
    std::stringstream ss;
    ss << pt.get<std::string>("Method.axis1");
    ss >> xis[0].x >> xis[0].y;
    xis[1] = xis[0];
    ss.str("");
    ss.clear();
    ss << pt.get<std::string>("Method.axis2");
    ss >> xis[2].x >> xis[2].y;
    xis[3] = xis[2];
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


