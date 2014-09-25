#pragma once
#include <math.h>
#include "vector2.h"


class tensor2d
{
    public:

        double   t00;
        double   t01;
        double   t10;
        double   t11;

        tensor2d() {}

        tensor2d(double p00, double p01, double p10, double p11)
        {
            t00 = p00;
            t01 = p01;
            t10 = p10;
            t11 = p11;
        }

        tensor2d& Set(double p00, double p01, double p10, double p11)
        {
            t00 = p00;
            t01 = p01;
            t10 = p10;
            t11 = p11;
            return (*this);
        }

        tensor2d& operator +=(const tensor2d& v)
        {
            t00 += v.t00;
            t01 += v.t01;
            t10 += v.t10;
            t11 += v.t11;
            return (*this);
        }

        tensor2d& operator -=(const tensor2d& v)
        {
            t00 -= v.t00;
            t01 -= v.t01;
            t10 -= v.t10;
            t11 -= v.t11;
            return (*this);
        }

        tensor2d& operator *=(double t)
        {
            t00 *= t;
            t01 *= t;
            t10 *= t;
            t11 *= t;
            return (*this);
        }

        tensor2d& operator /=(double t)
        {
            double f = 1.0F / t;
            t00 *= f;
            t01 *= f;
            t10 *= f;
            t11 *= f;
            return (*this);
        }

        tensor2d operator -(void) const
        {
            return (tensor2d(-t00, -t01, -t10, -t11));
        }

        tensor2d operator +(void) const
        {
            return (tensor2d(+t00, +t01, +t10, +t11));
        }

        tensor2d operator +(const tensor2d& v) const
        {
            return (tensor2d(t00 + v.t00, t01 + v.t01, t10 + v.t10, t11 + v.t11));
        }

        tensor2d operator -(const tensor2d& v) const
        {
            return (tensor2d(t00 - v.t00, t01 - v.t01, t10 - v.t10, t11 - v.t11));
        }

        tensor2d operator *(double t) const
        {
            return (tensor2d(t00 * t, t01 * t, t10 * t, t11 * t));
        }

        tensor2d operator /(double t) const
        {
            double f = 1.0F / t;
            return (tensor2d(t00 * f, t01 * f, t10 * f, t11 * f));
        }

        double operator *(const tensor2d& v) const
        {
            return (t00 * v.t00 + t01 * v.t01 + t10 * v.t10 + t11 * v.t11);
        }

        vector2d operator *(const vector2d& v) const
        {
            return vector2d(t00 * v.x + t01 * v.y, t10 * v.x + t11 * v.y);
        }

};

tensor2d vector2d::operator ^(const vector2d& v) const
{
    return (tensor2d(x * v.x, x * v.y, y * v.x, y * v.y));
}

inline tensor2d operator *(double t, const tensor2d& T)
{
    return (tensor2d(t * T.t00, t * T.t01, t * T.t10, t * T.t11));
}
