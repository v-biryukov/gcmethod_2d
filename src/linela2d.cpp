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
#include <algorithm>
#include <sstream>
#include <omp.h>

//#define USE_SIMPLE_OMEGA

// ////////////////////////////////////////////////////////////
#ifdef USE_SIMPLE_OMEGA
// /////////////////////////////////////////////////////////

riemann_data linela2d::get_riemann_inv_X(int pn, int en0, int en1)
{
    point_data pd = data[pn];
    riemann_data rd;


    rd.w[0] =  0;//-c30/c10 * pd.sxx   +     pd.syy;

    rd.w[1] = -0.5 * c2[en0]*rho[en0]  * pd.vy    +   0.5                     * pd.sxy;

    rd.w[2] = +0.5 * c2[en1]*rho[en1]  * pd.vy    +   0.5                     * pd.sxy;

    rd.w[3] = - 0.5 * rho[en0] * c3(en0) * pd.vx    +   0.5 * c3(en0) / c1[en0] * pd.sxx;

    rd.w[4] = + 0.5 * rho[en1] * c3(en1) * pd.vx    +   0.5 * c3(en1) / c1[en1] * pd.sxx;

    return rd;
}


riemann_data linela2d::get_riemann_inv_X(int pn, int en)
{
    return get_riemann_inv_X(pn, en, en);
}

void linela2d::set_point_data_X(int pn, riemann_data & rd)
{
    int en1 = eldata_X[pn].el[0] > 0 ? eldata_X[pn].el[0] : 0;
    int en2 = eldata_X[pn].el[1] > 0 ? eldata_X[pn].el[1] : 0;
    point_data pd;

    pd.vx = - 1.0/(rho[en1] * c3(en1)) * rd.w[3] + 1.0/(rho[en2] * c3(en2)) * rd.w[4];

    pd.vy = - 1.0 / (rho[en1] * c2[en1]) * rd.w[1] + 1.0 / (rho[en2] * c2[en2]) * rd.w[2];

    pd.sxx = +c1[en1]/c3(en1) * rd.w[3] + c1[en2]/c3(en2) * rd.w[4];

    pd.sxy = rd.w[1] + rd.w[2];

    pd.syy = rd.w[0] + rd.w[3] + rd.w[4];

    data_new[pn] = data[pn] + pd;
}



riemann_data linela2d::get_riemann_inv_Y(int pn, int en0, int en1)
{
    point_data pd = data[pn];
    riemann_data rd;

    rd.w[0] =  0.0;//1.0                    * pd.sxx   - c30 / c10 * pd.syy;

    rd.w[1] = -0.5 * c2[en0]*rho[en0] * pd.vx    +   0.5     * pd.sxy;

    rd.w[2] = +0.5 * c2[en1]*rho[en1] * pd.vx    +   0.5     * pd.sxy;

    rd.w[3] = -0.5 * c1[en0]*rho[en0] * pd.vy    +   0.5     * pd.syy;

    rd.w[4] = +0.5 * c1[en1]*rho[en1] * pd.vy    +   0.5     * pd.syy;
    return rd;
}


riemann_data linela2d::get_riemann_inv_Y(int pn, int en)
{
    return get_riemann_inv_Y(pn, en, en);
}

void linela2d::set_point_data_Y(int pn, riemann_data & rd)
{
    int en1 = eldata_Y[pn].el[0] > 0 ? eldata_Y[pn].el[0] : 0;
    int en2 = eldata_Y[pn].el[1] > 0 ? eldata_Y[pn].el[1] : 0;
    point_data pd;

    pd.vx = - 1.0 / (c2[en1] * rho[en1]) * rd.w[1] + 1.0 / (c2[en2] * rho[en2]) * rd.w[2];

    pd.vy = - 1.0 / (c1[en1] * rho[en1]) * rd.w[3] + 1.0 / (c1[en2] * rho[en2]) * rd.w[4];

    pd.sxx = rd.w[0] + c3(en1)/c1[en1] * rd.w[3] + c3(en2)/c1[en2] * rd.w[4];

    pd.sxy = rd.w[1] + rd.w[2];

    pd.syy = rd.w[3] + rd.w[4];

    data_new[pn] = data[pn] + pd;
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



// ////////////////////////////////////////////////////////////
#endif
// /////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////
#ifndef USE_SIMPLE_OMEGA
// /////////////////////////////////////////////////////////



riemann_data linela2d::get_riemann_inv_X(int pn, int en0, int en1)
{
    point_data pd = data[pn];
    riemann_data rd;

    double c1R = c1[en0];
    double c1L = c1[en1];
    double c2R = c2[en0];
    double c2L = c2[en1];
    double rR = rho[en0];
    double rL = rho[en1];
    double c3R = c3(en0);
    double c3L = c3(en1);

    double c30 = (c3L+c3R)/2.0;
    double c10 = (c1L+c1R)/2.0;

    rd.w[0] =  -c30/c10 * pd.sxx + pd.syy;

    rd.w[1] =  -c3R*rR/2.0 * pd.vx + c3R/c1R/2.0*pd.sxx;

    rd.w[2] =  +c3L*rL/2.0 * pd.vx + c3L/c1L/2.0*pd.sxx;

    rd.w[3] =  -c2R*rR/2.0*pd.vy + 0.5 * pd.sxy;

    rd.w[4] =  +c2L*rL/2.0*pd.vy + 0.5 * pd.sxy;

    return rd;
}

riemann_data linela2d::get_riemann_inv_X(int pn, int en)
{
    return get_riemann_inv_X(pn, en, en);
}

void linela2d::set_point_data_X(int pn, riemann_data & rd)
{
    int en0 = eldata_X[pn].el[0] > 0 ? eldata_X[pn].el[0] : 0;
    int en1 = eldata_X[pn].el[1] > 0 ? eldata_X[pn].el[1] : 0;
    point_data pd;

    double c1R = c1[en0];
    double c1L = c1[en1];
    double c2R = c2[en0];
    double c2L = c2[en1];
    double rR = rho[en0];
    double rL = rho[en1];
    double c3R = c3(en0);
    double c3L = c3(en1);

    double c30 = (c3L+c3R)/2.0;
    double c10 = (c1L+c1R)/2.0;

    pd.vx =  2.0*(-c1R/c3R * rd.w[1] + c1L/c3L * rd.w[2])/(c1L*rL+c1R*rR);

    pd.vy =  2.0*(-rd.w[3] + rd.w[4])/(c2L*rL+c2R*rR);

    pd.sxx = 2.0*c1L*c1R*(rL/c3R * rd.w[1] + rR/c3L * rd.w[2])/(c1L*rL+c1R*rR);

    pd.sxy = 2.0*(c2L*rL * rd.w[3] + c2R*rR * rd.w[4])/(c2L*rL+c2R*rR);

    pd.syy = rd.w[0] + 2.0*(rL/c3R * rd.w[1] + rR/c3L * rd.w[2])*c1L*c1R*c30/(c1L*rL+c1R*rR)/c10;

    data_new[pn] = data[pn] + pd;
}

riemann_data linela2d::get_riemann_inv_Y(int pn, int en0, int en1)
{
    point_data pd = data[pn];
    riemann_data rd;

    double c1R = c1[en0];
    double c1L = c1[en1];
    double c2R = c2[en0];
    double c2L = c2[en1];
    double rR = rho[en0];
    double rL = rho[en1];
    double c3R = c3(en0);
    double c3L = c3(en1);

    double c30 = (c3L+c3R)/2.0;
    double c10 = (c1L+c1R)/2.0;


    rd.w[0] = pd.sxx - c30/c10 * pd.syy;

    rd.w[1] = -c1R*rR/2.0 * pd.vy + 0.5 * pd.syy;

    rd.w[2] = +c1L*rL/2.0 * pd.vy + 0.5 * pd.syy;

    rd.w[3] = -c2R*rR/2.0 * pd.vx + 0.5 * pd.sxy;

    rd.w[4] = +c2L*rL/2.0 * pd.vx + 0.5 * pd.sxy;

    return rd;
}
riemann_data linela2d::get_riemann_inv_Y(int pn, int en)
{
    return get_riemann_inv_Y(pn, en, en);
}

void linela2d::set_point_data_Y(int pn, riemann_data & rd)
{
    int en0 = eldata_Y[pn].el[0] > 0 ? eldata_Y[pn].el[0] : 0;
    int en1 = eldata_Y[pn].el[1] > 0 ? eldata_Y[pn].el[1] : 0;

    point_data pd;

    double c1R = c1[en0];
    double c1L = c1[en1];
    double c2R = c2[en0];
    double c2L = c2[en1];
    double rR = rho[en0];
    double rL = rho[en1];
    double c3R = c3(en0);
    double c3L = c3(en1);

    double c30 = (c3L+c3R)/2.0;
    double c10 = (c1L+c1R)/2.0;


    pd.vx = 2.0*(-rd.w[3] + rd.w[4])/(c2L*rL+c2R*rR);

    pd.vy = 2.0*(-rd.w[1] + rd.w[2])/(c1L*rL+c1R*rR);

    pd.sxx = rd.w[0] + 2.0*c30/c10*(c1L*rL * rd.w[1] + c1R*rR * rd.w[2])/(c1L*rL+c1R*rR);

    pd.sxy = 2.0*(c2L*rL * rd.w[3] + c2R*rR * rd.w[4])/(c2L*rL+c2R*rR);

    pd.syy = 2.0*(c1L*rL * rd.w[1] + c1R*rR * rd.w[2])/(c1L*rL+c1R*rR);

    data_new[pn] = data[pn] + pd;
}

std::vector<double> linela2d::get_lambda_X(int point_n)
{
    std::vector<double> lambda(5);
    lambda[0] = 0.0;
    lambda[1] = - c1[eldata_X[point_n].el[0]];
    lambda[2] = + c1[eldata_X[point_n].el[1]];
    lambda[3] = - c2[eldata_X[point_n].el[0]];
    lambda[4] = + c2[eldata_X[point_n].el[1]];
    return lambda;
}

std::vector<double> linela2d::get_lambda_Y(int point_n)
{
    std::vector<double> lambda(5);
    lambda[0] = 0.0;
    lambda[1] = - c1[eldata_Y[point_n].el[0]];
    lambda[2] = + c1[eldata_Y[point_n].el[1]];
    lambda[3] = - c2[eldata_Y[point_n].el[0]];
    lambda[4] = + c2[eldata_Y[point_n].el[1]];
    return lambda;
}


// ////////////////////////////////////////////////////////////
#endif
// /////////////////////////////////////////////////////////

void linela2d::calculate_point_elements()
{
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        for (int j = 0; j < 2; j++)
        {
            eldata_X[i].el[j] = -1;
            eldata_Y[i].el[j] = -1;
        }
    }
    for (int i = 0; i < mesh->get_number_of_points(); i++)
        for (int k = 1; k < 5; k++)
        {
            vector2d p = mesh->get_point(i) - lambda_hint * tau * directions[k];
            if (mesh->point_types[i] == mesh->INNER || mesh->is_inside_contour(p, i))
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
                        break;
                    }
                }
        }
}


// TODO
void step_any(vector2d dir)
{
}


void linela2d::step_X()
{
    #pragma omp parallel for
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        std::vector<std::vector<riemann_data> > rdata;
        rdata.resize(2, std::vector<riemann_data>((N+1)*(N+2)/2));
        int el0 = eldata_X[i].el[0];
        int el1 = eldata_X[i].el[1];
        for (int j = 0; j < (N+1)*(N+2)/2; j++)
        {
            if (el0 >= 0)
                rdata[0][j] = get_riemann_inv_X(mesh->elements[el0][j], el0);
            if (el1 >= 0)
                rdata[1][j] = get_riemann_inv_X(mesh->elements[el1][j], el1);
        }
        riemann_data rdata_here = get_riemann_inv_X(i, el0, el1);
        std::vector<double> lambda = get_lambda_X(i);
        riemann_data diff_rd = riemann_data();
        for (int k = 1; k < 5; k++)
        {
            vector2d p = mesh->get_point(i) - lambda[k] * tau * directions[2];
            int parity = (k-1)%2;
            if (eldata_X[i].el[parity] >= 0)
                diff_rd.w[k] = approximate(p, eldata_X[i].el[parity], rdata[parity] , k) - rdata_here.w[k];
        }

        set_point_data_X(i, diff_rd);
    }
    postprocess_border_conditions(0);
    data_new.swap(data);
}

void linela2d::step_Y()
{
    #pragma omp parallel for
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        std::vector<std::vector<riemann_data> > rdata;
        rdata.resize(2, std::vector<riemann_data>((N+1)*(N+2)/2));
        int el0 = eldata_Y[i].el[0];
        int el1 = eldata_Y[i].el[1];
        for (int j = 0; j < (N+1)*(N+2)/2; j++)
        {
            if (el0 >= 0)
                rdata[0][j] = get_riemann_inv_Y(mesh->elements[el0][j], el0);
            if (el1 >= 0)
                rdata[1][j] = get_riemann_inv_Y(mesh->elements[el1][j], el1);
        }
        riemann_data rdata_here = get_riemann_inv_Y(i, el0, el1);
        std::vector<double> lambda = get_lambda_Y(i);
        riemann_data diff_rd = riemann_data();
        for (int k = 1; k < 5; k++)
        {
            vector2d p = mesh->get_point(i) - lambda[k] * tau * directions[4];
            int parity = (k-1)%2;
            if (eldata_Y[i].el[parity] >= 0)
                diff_rd.w[k] = approximate(p, eldata_Y[i].el[parity], rdata[parity] , k) - rdata_here.w[k];
        }
        set_point_data_Y(i, diff_rd);
    }
    postprocess_border_conditions(1);
    data_new.swap(data);
}

void linela2d::postprocess_border_conditions(int axis)
{
    vector2d n0 = directions[2];
    vector2d n1 = directions[4];
    vector2d n;
    if (axis == 0)
        n = n0;
    else if (axis == 1)
        n = n1;
    else
    {
        std::cerr << "Error in postprocess_border_conditions" << std::endl;
        std::exit(1);
    }
    for (int i = 0; i < mesh->get_number_of_points(); i++)
    {
        if ((axis == 0 && (eldata_X[i].el[0] == -1 || eldata_X[i].el[1] == -1)) ||
            (axis == 1 && (eldata_Y[i].el[0] == -1 || eldata_Y[i].el[1] == -1)) )
        {
            if (mesh->point_types[i] == mesh_2d::ABSORB)
            {
                // do nothing
            }
            else if (mesh->point_types[i] == mesh_2d::FREE || mesh->point_types[i] == mesh_2d::FORCE)
            {
                double c1i=0.0, c2i=0.0, rhoi=0.0;
                for (int j = 0; j < mesh->triangles[i].size(); j++)
                {
                    c1i += c1[mesh->triangles[i][j]];
                    c2i += c2[mesh->triangles[i][j]];
                    rhoi += rho[mesh->triangles[i][j]];
                }
                c1i /= mesh->triangles[i].size();
                c2i /= mesh->triangles[i].size();

                double c3i = c1i - 2.0*c2i*c2i/c1i;
                rhoi /= mesh->triangles[i].size();

                vector2d p = mesh->point_normals[i];

                vector2d force;
                if (mesh->point_types[i] == mesh_2d::FORCE)
                {
                    double f = 40.0;
                    double A = 5.0;
                    double delay = 3.0*sqrt(1.5)/M_PI/f;
                    double X = M_PI*f*(time - delay);
                    force = - p * A *(1 - 2*X*X)*exp(-X*X);

                    //force = - p * A * sin( 2*M_PI * time / mesh->get_courant_time_step() / 53.0 );
                }
                else
                {
                    force = vector2d(0.0, 0.0);
                }

    //            // Simplified free border
    //            vector2d z = vector2d(data[i].sxx*p.x + data[i].sxy*p.y, data[i].sxy*p.x + data[i].syy*p.y) - force;
    //            double lambda = c1i*c1i*rhoi - 2.0*c2i*c2i*rhoi;
    //            double mu     = c2i*c2i*rhoi;
    //            point_data u;
    //            u.vx =  - 1.0/(rhoi*c2i)*z.x + 1.0/rhoi*(1.0/c2i-1.0/c1i)*(z*p)*p.x;
    //            u.vy =  - 1.0/(rhoi*c2i)*z.y + 1.0/rhoi*(1.0/c2i-1.0/c1i)*(z*p)*p.y;
    //            u.sxx = - 2.0*(z.x*p.x)      - (z*p)/(lambda + 2*mu) * (lambda - 2.0*(lambda+mu)*p.x*p.x);
    //            u.sxy = - (z.x*p.y+z.y*p.x)  - (z*p)/(lambda + 2*mu) * (        -2.0*(lambda+mu)*p.x*p.y);
    //            u.syy = - 2.0*(z.y*p.y)      - (z*p)/(lambda + 2*mu) * (lambda - 2.0*(lambda+mu)*p.y*p.y);
    //            data_new[i] += u;

                vector2d z = vector2d(data[i].sxx*p.x + data[i].sxy*p.y, data[i].sxy*p.x + data[i].syy*p.y) - force;
                double om1 = (2.0*(p*n)*(n*z) - (p*z))/((c1i+c3i)*(n*p)*(n*p)-c3i*(p*p));
                vector2d b = (z-om1*c3i*p)/(c2i*(n*p));



                double eps = 1e-6;
                point_data u = {0.0, 0.0, 0.0, 0.0, 0.0};
                if (n*p > eps)
                {
                    u.vx =  +(om1*n.x - (n*b)*n.x + b.x)/rhoi;
                    u.vy =  +(om1*n.y - (n*b)*n.y + b.y)/rhoi;
                }
                else if (n*p < -eps)
                {
                    u.vx =  -(om1*n.x - (n*b)*n.x + b.x)/rhoi;
                    u.vy =  -(om1*n.y - (n*b)*n.y + b.y)/rhoi;
                }

                u.sxx = -((c1i-c3i)*om1-2*c2i*(n*b))*n0.x*n0.x - c3i*om1 - 2.0*c2i*(n.x*b.x);
                u.sxy = -((c1i-c3i)*om1-2*c2i*(n*b))*n0.x*n0.y           - c2i*(n.x*b.y+n.y*b.x);
                u.syy = -((c1i-c3i)*om1-2*c2i*(n*b))*n0.y*n0.y - c3i*om1 - 2.0*c2i*(n.y*b.y);

                data_new[i] += u;
            }
        }
    }
}

void linela2d::step()
{
    step_X();
    step_Y();
    time += tau;
}

void linela2d::calculate()
{
    print_parameters();
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

double linela2d::approximate(vector2d r, int tn, std::vector<riemann_data> & rdata, int k)
{
    if (!mesh->is_inside(r, tn))
    {
        //std::cerr << "Error, point not in the cell number " << tn << "!\n";
        //std::exit(1);
    }
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
        v += sa * rdata[0].w[k];
        v += sb * rdata[1].w[k];
        v += sc * rdata[2].w[k];

    }
    else if ( N == 2 )
    {
        v += sa*(2*sa-1) * rdata[0].w[k];
        v += sb*(2*sb-1) * rdata[1].w[k];
        v += sc*(2*sc-1) * rdata[2].w[k];
        v += 4*sa*sb     * rdata[3].w[k];
        v += 4*sb*sc     * rdata[4].w[k];
        v += 4*sc*sa     * rdata[5].w[k];
    }
    else if ( N == 3 )
    {
        v += sa*(3*sa-1)*(3*sa-2)/2 * rdata[0].w[k];
        v += sb*(3*sb-1)*(3*sb-2)/2 * rdata[1].w[k];
        v += sc*(3*sc-1)*(3*sc-2)/2 * rdata[2].w[k];
        v += 9*sa*(3*sa-1)*sb/2     * rdata[3].w[k];
        v += 9*sb*(3*sb-1)*sa/2     * rdata[4].w[k];
        v += 9*sb*(3*sb-1)*sc/2     * rdata[5].w[k];
        v += 9*sc*(3*sc-1)*sb/2     * rdata[6].w[k];
        v += 9*sc*(3*sc-1)*sa/2     * rdata[7].w[k];
        v += 9*sa*(3*sa-1)*sc/2     * rdata[8].w[k];
        v += 27*sa*sb*sc            * rdata[9].w[k];
    }
    else if ( N == 4 )
    {
        v += sa*(4*sa-1)*(2*sa-1)*(4*sa-3)/3 * rdata[0].w[k];
        v += sb*(4*sb-1)*(2*sb-1)*(4*sb-3)/3 * rdata[1].w[k];
        v += sc*(4*sc-1)*(2*sc-1)*(4*sc-3)/3 * rdata[2].w[k];
        v += 16*sa*(4*sa-1)*(2*sa-1)*sb/3    * rdata[3].w[k];
        v += 4*sa*(4*sa-1)*sb*(4*sb-1)       * rdata[4].w[k];
        v += 16*sb*(4*sb-1)*(2*sb-1)*sa/3    * rdata[5].w[k];
        v += 16*sb*(4*sb-1)*(2*sb-1)*sc/3    * rdata[6].w[k];
        v += 4*sb*(4*sb-1)*sc*(4*sc-1)       * rdata[7].w[k];
        v += 16*sc*(4*sc-1)*(2*sc-1)*sb/3    * rdata[8].w[k];
        v += 16*sc*(4*sc-1)*(2*sc-1)*sa/3    * rdata[9].w[k];
        v += 4*sa*(4*sa-1)*sc*(4*sc-1)       * rdata[10].w[k];
        v += 16*sa*(4*sa-1)*(2*sa-1)*sc/3    * rdata[11].w[k];
        v += 32*sa*(4*sa-1)*sb*sc            * rdata[12].w[k];
        v += 32*sb*(4*sb-1)*sc*sa            * rdata[13].w[k];
        v += 32*sc*(4*sc-1)*sa*sb            * rdata[14].w[k];
    }
    return v;
}

void linela2d::set_directions(double angle)
{
    directions.resize(5);
    directions[0] = vector2d(0.0, 0.0);
    directions[2] = vector2d(cos(angle), sin(angle));
    directions[1] = - directions[2];
    directions[4] = vector2d(-sin(angle), cos(angle));
    directions[3] = - directions[4];
}


void linela2d::print_parameters()
{
    using std::cout;
    using std::endl;

    cout << "Start linear elastics calculation using grid-characteristic method" << std::endl;
    cout << "Parameters:" << std::endl;
    cout << "\t Number of main points:     " << mesh->get_number_of_main_points() << endl;
    cout << "\t Number of all points:      " << mesh->get_number_of_points() << endl;
    cout << "\t Number of elements:        " << mesh->get_number_of_elements() << endl;
    cout << "\t Number of all triangles:   " << mesh->get_number_of_elements() * N * N << endl;
    cout << "\t Space step (approximately):" << mesh->get_h() << endl;
    cout << "\t Min element altitude:      " << mesh->get_min_altitude() << endl;
    cout << "\t Min triangle altitude:     " << mesh->get_min_altitude()/N << endl;
    cout << "\t Courant time step:         " << mesh->get_courant_time_step() << endl;
    cout << "\t Courant number:            " << courant_multiplier << endl;
    cout << "\t Time step:                 " << tau << endl;
    if (courant_multiplier > 1.0)
    {
        cout << endl;
        cout << "\t\t WARNING! Courant number > 1";
        cout << endl;
    }
    cout << "\t Number of time steps:      " << number_of_steps << endl;
    cout << "\t Order of approximation:    " << N;
    if (is_monotonic)
        cout << ", monotonic";
    cout << endl;
    cout << "\t Longitudinal speed of sound vary from " << *std::min_element(c1.begin(), c1.end())   << " to " << *std::max_element(c1.begin(), c1.end()) << endl;
    cout << "\t Transversal speed of sound vary from  " << *std::min_element(c2.begin(), c2.end())   << " to " << *std::max_element(c2.begin(), c2.end()) << endl;
    cout << "\t Density vary from                     " << *std::min_element(rho.begin(), rho.end()) << " to " << *std::max_element(rho.begin(), rho.end()) << endl;
    cout << endl;
}

void linela2d::init()
{
    time = 0.0;
    c1.resize(mesh->get_number_of_elements());
    c2.resize(mesh->get_number_of_elements());
    rho.resize(mesh->get_number_of_elements());
    data.resize(mesh->get_number_of_points());
    data_new.resize(mesh->get_number_of_points());
    eldata_X.resize(mesh->get_number_of_points());
    eldata_Y.resize(mesh->get_number_of_points());
    lambda_hint = mesh->get_min_sound_speed();

    set_directions(0.0);

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


    string time_step_choosing_type = pt.get<string>("Method.time_step_choosing_type");
    if (time_step_choosing_type == "courant")
    {
        courant_multiplier = pt.get<double>("Method.courant_multiplier");
        tau = mesh->get_courant_time_step() * courant_multiplier;
    }
    else if (time_step_choosing_type == "tau")
    {
        tau = pt.get<double>("Method.tau");
        courant_multiplier = tau / mesh->get_courant_time_step();
    }
}

void linela2d::save_to_vtk(std::string name)
{
    std::ofstream vtk_file(name.c_str());
    vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
    vtk_file << "DATASET POLYDATA\nPOINTS " << mesh->get_number_of_points() << " float\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
        vtk_file << mesh->get_point(i).x << " " << mesh->get_point(i).y << " "  << 0.0 << "\n";
    vtk_file << "\nPOLYGONS " << mesh->get_number_of_elements() << " " << mesh->get_number_of_elements()*4 << "\n";
    for (int i = 0; i < mesh->get_number_of_elements(); i++)
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
    vtk_file << "SCALARS rho FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        double result = 0.0;
        if (!mesh->triangles[i].empty())
        {
            for (int j = 0; j < mesh->triangles[i].size(); j++)
                result += rho[mesh->triangles[i][j]];
            result /= mesh->triangles[i].size();
        }
        vtk_file <<  result << "\n";

    }
    vtk_file << "SCALARS c1 FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        double result = 0.0;
        if (!mesh->triangles[i].empty())
        {
            for (int j = 0; j < mesh->triangles[i].size(); j++)
                result += c1[mesh->triangles[i][j]];
            result /= mesh->triangles[i].size();
        }
        vtk_file <<  result << "\n";
    }
    vtk_file << "SCALARS c2 FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( int i = 0; i < mesh->get_number_of_points(); i++ )
    {
        double result = 0.0;
        if (!mesh->triangles[i].empty())
        {
            for (int j = 0; j < mesh->triangles[i].size(); j++)
                result += c2[mesh->triangles[i][j]];
            result /= mesh->triangles[i].size();
        }
        vtk_file <<  result << "\n";
    }
}
