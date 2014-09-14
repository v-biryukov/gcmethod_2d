#include "gcmethod_2d.h"
#include <math.h>


int main()
{
    mesh_2d m = mesh_2d("gcmethod_2d.ini");
    m.create_mesh();
    gcmethod_2d g = gcmethod_2d("gcmethod_2d.ini", m);
    g.calculate();
    std::ofstream f ("results.txt");
    f << std::fixed << std::setprecision(12) << g.L(1) << '\n' << g.L(2) << '\n' << g.L_inf();;
}

/*
int main()
{
    std::vector<double> hs = { 0.8, 0.4, 0.2, 0.1, 0.05 };
    std::vector<std::vector<double> > L1s, L2s, Linfs;
    L1s.resize(hs.size());
    L2s.resize(hs.size());
    Linfs.resize(hs.size());
    std::ofstream f ("out/errors.txt");
	for ( int i = 0; i < hs.size(); i++ )
	{
        mesh_2d m = mesh_2d("gcmethod_2d.ini");
        m.set_h(hs.at(i));
        std::cout << "\nCreating mesh with h = " << m.get_h() << "\n";

        m.create_mesh();
        for ( int N = 1; N <= 4; N++ )
        {
            gcmethod_2d test = gcmethod_2d("gcmethod_2d.ini", m);
            std::cout << "\nCalculating task with h = " << m.get_h() << ", N = " << N << "\n";

            test.N = N;
            test.calculate();
            L1s.at(i).push_back(test.L(1));
            L2s.at(i).push_back(test.L(2));
            Linfs.at(i).push_back(test.L_inf());
        }
    }
    for ( int N = 1; N <= 4; N++ )
    {
        for ( int i = 0; i < hs.size(); i++ )
        {
            f << "\n\n  Task with h = " << std::fixed << std::setprecision(3) << hs.at(i) << ", N = " << N << " :\n";
            f << "    L1 = "<< std::fixed << std::setprecision(12) << L1s.at(i).at(N-1) << " \n";
            f << "    L2 = "<< std::fixed << std::setprecision(12) << L2s.at(i).at(N-1) << " \n";
            f << "    Linf = "<< std::fixed << std::setprecision(12) << Linfs.at(i).at(N-1) << " \n";
            if ( i > 0 )
            {
                f << "    P_1 = log(L1(h) - L1(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(L1s.at(i).at(N-1)) - log(L1s.at(i-1).at(N-1)))/log(2.0) << "\n";
                f << "    P_2 = log(L2(h) - L2(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(L2s.at(i).at(N-1)) - log(L2s.at(i-1).at(N-1)))/log(2.0) << "\n";
                f << "    P_inf = log(Linf(h) - Linf(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(Linfs.at(i).at(N-1)) - log(Linfs.at(i-1).at(N-1)))/log(2.0) << "\n";
            }
        }

    }
	return 0;
}
*/
