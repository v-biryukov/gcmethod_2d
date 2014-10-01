#pragma once
#include <math.h>


class tensor2d;
class vector2d
{
    public:

        double   x;
        double   y;

        vector2d() {}

        vector2d(double r, double s)
        {
            x = r;
            y = s;
        }

        vector2d& Set(double r, double s)
        {
            x = r;
            y = s;
            return (*this);
        }

        double& operator [](long k)
        {
            return ((&x)[k]);
        }

        const double& operator [](long k) const
        {
            return ((&x)[k]);
        }

        vector2d& operator +=(const vector2d& v)
        {
            x += v.x;
            y += v.y;
            return (*this);
        }

        vector2d& operator -=(const vector2d& v)
        {
            x -= v.x;
            y -= v.y;
            return (*this);
        }

        vector2d& operator *=(double t)
        {
            x *= t;
            y *= t;
            return (*this);
        }

        vector2d& operator /=(double t)
        {
            double f = 1.0F / t;
            x *= f;
            y *= f;
            return (*this);
        }

        vector2d& operator &=(const vector2d& v)
        {
            x *= v.x;
            y *= v.y;
            return (*this);
        }

        vector2d operator -(void) const
        {
            return (vector2d(-x, -y));
        }

        vector2d operator +(const vector2d& v) const
        {
            return (vector2d(x + v.x, y + v.y));
        }

        vector2d operator -(const vector2d& v) const
        {
            return (vector2d(x - v.x, y - v.y));
        }

        vector2d operator *(double t) const
        {
            return (vector2d(x * t, y * t));
        }

        vector2d operator /(double t) const
        {
            double f = 1.0F / t;
            return (vector2d(x * f, y * f));
        }

        double operator *(const vector2d& v) const
        {
            return (x * v.x + y * v.y);
        }

        vector2d operator &(const vector2d& v) const
        {
            return (vector2d(x * v.x, y * v.y));
        }

        bool operator ==(const vector2d& v) const
        {
            return ((fabs(x - v.x) < 1e-10) && (fabs(y - v.y) < 1e-10));
        }

        bool operator !=(const vector2d& v) const
        {
            return ((fabs(x - v.x) > 1e-10) || (fabs(y - v.y) > 1e-10));
        }

        vector2d& Normalize(void)
        {
            return (*this /= sqrtf(x * x + y * y));
        }

		double angle(void) const
        {
            return atan2(this->y, this->x);
        }
        tensor2d operator ^(const vector2d& v) const;
};



class point2d : public vector2d
{
    public:

        point2d() {}

        point2d(double r, double s) : vector2d(r, s) {}

        point2d& operator =(const vector2d& v)
        {
            x = v.x;
            y = v.y;
            return (*this);
        }

        point2d& operator *=(double t)
        {
            x *= t;
            y *= t;
            return (*this);
        }

        point2d& operator /=(double t)
        {
            double f = 1.0F / t;
            x *= f;
            y *= f;
            return (*this);
        }

        point2d operator -(void) const
        {
            return (point2d(-x, -y));
        }

        point2d operator +(const vector2d& v) const
        {
            return (point2d(x + v.x, y + v.y));
        }

        point2d operator -(const vector2d& v) const
        {
            return (point2d(x - v.x, y - v.y));
        }

        vector2d operator -(const point2d& p) const
        {
            return (vector2d(x - p.x, y - p.y));
        }

        point2d operator *(double t) const
        {
            return (point2d(x * t, y * t));
        }

        point2d operator /(double t) const
        {
            double f = 1.0F / t;
            return (point2d(x * f, y * f));
        }
};


inline vector2d operator *(double t, const vector2d& v)
{
    return (vector2d(t * v.x, t * v.y));
}

inline point2d operator *(double t, const point2d& p)
{
    return (point2d(t * p.x, t * p.y));
}

inline double Dot(const vector2d& v1, const vector2d& v2)
{
    return (v1 * v2);
}

inline double Vec(const vector2d& v1, const vector2d& v2)
{
    return (v1.x * v2.y - v2.x * v1.y);
}

inline double Magnitude(const vector2d& v)
{
    return (sqrt(v.x * v.x + v.y * v.y));
}

inline double InverseMag(const vector2d& v)
{
    return (1.0F / sqrtf(v.x * v.x + v.y * v.y));
}

inline double SquaredMag(const vector2d& v)
{
    return (v.x * v.x + v.y * v.y);
}
