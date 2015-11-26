#include <cmath>

//3-Dimensional vector class
template<typename T>
class Vec3
{
public:
    //The 3D coordinates
    T x, y, z;

    //Empty constructor
    Vec3()
    {
        x=0;
        y=0;
        z=0;
    }

    //One coordinate constructor
    Vec3(T xx)
    {
        x=xx;
        y=xx;
        z=xx;
    }

    //Full constructor
    Vec3(T xx, T yy, T zz)
    {
        x=xx;
        y=yy;
        z=zz;
    }

    //Calculate length (magnitude) of the vector
    T length()
    {
        return sqrt(x*x + y*y + z*z);
    }

    //Normalize the vector
    Vec3<T>& normalize()
    {
        T len = length();
        if( len>0 )
        {
            T invLen = 1/len;

            x*=invLen;
            y*=invLen;
            z*=invLen;
        }

        return *this;
    }

    //Dot product of the vector
    T dot(Vec3<T> v)
    {
        return x*v.x + y*v.y + z*v.z;
    }

    //Cross product of the vector
    Vec3<T> cross(Vec3<T> v)
    {
        return Vec3<T>(y*v.z - z*v.y,
                       z*v.x - x*v.z,
                       x*v.y - y*v.x);
    }

    //Vector addition
    Vec3<T> operator+( const Vec3<T> &v ) const
    {
        return Vec3<T>(x + v.x,
                       y + v.y,
                       z + v.z);
    }

    //Vector subtraction
    Vec3<T> operator-( const Vec3<T> &v ) const
    {
        return Vec3<T>(x - v.x,
                       y - v.y,
                       z - v.z);
    }

    //Scalar multiplication
    Vec3<T> operator*( const T &r ) const
    {
        return Vec3<T>(x*r,
                       y*r,
                       z*r);
    }
};
