#include "gcmethod_2d.h"
#include "mesh_2d.h"

class analyzer
{

    std::string results_path, ini_path;
    std::vector<std::vector<std::vector<double> > > L;
    std::vector<std::vector<std::vector<double> > > Lm;

public:
    analyzer(std::string ini_path_t, std::string results_path_t) : ini_path(ini_path_t), results_path(results_path_t) {};

    void analyze(std::vector<double> hs, std::vector<int> Ns = {1, 2, 3, 4})
    {
        L.resize(hs.size());
        Lm.resize(hs.size());
        std::ofstream f (results_path);
        for ( int i = 0; i < hs.size(); i++ )
        {
            mesh_2d m = mesh_2d("gcmethod_2d.ini");
            m.set_h(hs.at(i));
            std::cout << "\nCreating mesh with h = " << m.get_h() << "is_structured = " << m.get_is_structured() << "\n";
            m.create_mesh();
            for ( int N : Ns )
            {
                for ( bool is_monotonic : {true, false} )
                {
                    gcmethod_2d test = gcmethod_2d("gcmethod_2d.ini", m);
                    std::cout << "\nCalculating task with h = " << m.get_h() << ", N = " << N << ", is_monotonic = " << is_monotonic << "\n";
                    test.N = N;
                    test.is_monotonic = is_monotonic;
                    test.calculate();
                    if ( is_monotonic )
                        Lm.at(i).push_back({test.L(1), test.L(2), test.L_inf()});
                    else
                        L.at(i).push_back({test.L(1), test.L(2), test.L_inf()});
                }
            }
        }
        int n = 0;
        f << "   Non Monotonic   \n\n";
        for ( int N : Ns )
        {
            for ( int i = 0; i < hs.size(); i++ )
            {
                f << "\n\n  Task with h = " << std::fixed << std::setprecision(3) << hs.at(i) << ", N = " << N << " :\n";
                f << "    L1 = "<< std::fixed << std::setprecision(12) << L.at(i).at(n).at(0) << " \n";
                f << "    L2 = "<< std::fixed << std::setprecision(12) << L.at(i).at(n).at(1) << " \n";
                f << "    Linf = "<< std::fixed << std::setprecision(12) << L.at(i).at(n).at(2) << " \n";
                if ( i > 0 )
                {
                    f << "    P_1 = log(L1(h) - L1(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(L.at(i-1).at(n).at(0)) - log(L.at(i).at(n).at(0)))/log(2.0) << "\n";
                    f << "    P_2 = log(L2(h) - L2(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(L.at(i-1).at(n).at(1)) - log(L.at(i).at(n).at(1)))/log(2.0) << "\n";
                    f << "    P_inf = log(Linf(h) - Linf(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(L.at(i-1).at(n).at(2)) - log(L.at(i).at(n).at(2)))/log(2.0) << "\n";
                }
            }
            n++;
        }
        f << "\n\n   Monotonic   \n\n";
        n = 0;
        for ( int N : Ns )
        {
            for ( int i = 0; i < hs.size(); i++ )
            {
                f << "\n\n  Task with h = " << std::fixed << std::setprecision(3) << hs.at(i) << ", N = " << N << " :\n";
                f << "    L1 = "<< std::fixed << std::setprecision(12) << Lm.at(i).at(n).at(0) << " \n";
                f << "    L2 = "<< std::fixed << std::setprecision(12) << Lm.at(i).at(n).at(1) << " \n";
                f << "    Linf = "<< std::fixed << std::setprecision(12) << Lm.at(i).at(n).at(2) << " \n";
                if ( i > 0 )
                {
                    f << "    P_1 = log(L1(h) - L1(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(Lm.at(i-1).at(n).at(0)) - log(Lm.at(i).at(n).at(0)))/log(2.0) << "\n";
                    f << "    P_2 = log(L2(h) - L2(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(Lm.at(i-1).at(n).at(1)) - log(Lm.at(i).at(n).at(1)))/log(2.0) << "\n";
                    f << "    P_inf = log(Linf(h) - Linf(h/2))/log2 = " << std::fixed << std::setprecision(6) << -(log(Lm.at(i-1).at(n).at(2)) - log(Lm.at(i).at(n).at(2)))/log(2.0) << "\n";
                }
            }
            n++;
        }
    }


};

