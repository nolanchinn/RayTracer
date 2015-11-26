#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>

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
    T length() const
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

    //Access the vector coordinates as an array constant
    const T& operator[]( uint8_t i ) const
    {
        return (&x)[i];
    }

    //Access the vector coordinates as an array not constant
    T& operator[]( uint8_t i )
    {
        return (&x)[i];
    }

    //Vector output
    friend std::ostream& operator<<( std::ostream &os, const Vec3<T> &v )
    {
        return os << "<" << v.x << ", " << v.y << ", " << v.z << ">";
    }
};

//2 specialized versions, shorthand
typedef Vec3<float> Vec3f;
typedef Vec3<int> Vec3i;

//4x4 matrix class
template<typename T>
class Matrix44
{
public:
    T m[4][4] = {{1,0,0,0}, //the default matrix is the 4x4 identity matrix
                 {0,1,0,0},
                 {0,0,1,0},
                 {0,0,0,1}};

    //empty constructor, does nothing, but won't throw an error
    Matrix44(){}

    //constructor that defines all matrix positions
    Matrix44( T a, T b, T c, T d, T e, T f, T g, T h, T i, T j, T k, T l, T m, T n, T o, T p)
    {
        m[0][0] = a;
        m[0][1] = b;
        m[0][2] = c;
        m[0][3] = d;
        m[1][0] = e;
        m[1][1] = f;
        m[1][2] = g;
        m[1][3] = h;
        m[2][0] = i;
        m[2][1] = j;
        m[2][2] = k;
        m[2][3] = l;
        m[3][0] = m;
        m[3][1] = n;
        m[3][2] = o;
        m[3][3] = p;
    }

    //mat[i] instead of mat.m[i] (constant)
    const T* operator[]( uint8_t i ) const
    {
        return m[i];
    }

    //mat[i] instead of mat.m[i]
    T* operator[]( uint8_t i )
    {
        return m[i];
    }

    //matrix multiplication
    Matrix44 operator*( const Matrix44& rhs ) const
    {
        Matrix44 mult;
        for( uint8_t i = 0; i < 4; i++ )
        {
            for( uint8_t j = 0; j < 4; j++ )
            {
                mult[i][j] = m[i][0] * rhs[0][j] +
                             m[i][1] * rhs[1][j] +
                             m[i][2] * rhs[2][j] +
                             m[i][3] * rhs[3][j];
            }
        }

        return mult;
    }

    //matrix multiplication (c = a*b)
    static void multiply( const Matrix44<T> &a, const Matrix44& b, Matrix44 &c)
    {
        for( uint8_t i = 0; i < 4; i++ )
        {
            for( uint8_t j = 0; j < 4; j++ )
            {
                c[i][j] = a[i][0] * b[0][j] +
                          a[i][1] * b[1][j] +
                          a[i][2] * b[2][j] +
                          a[i][3] * b[3][j];
            }
        }
    }

    //transform a point using the matrix
    template<typename S>
    void multVecMatrix(const Vec3<S> &src, Vec3<S> &dst) const
    {
        S w;
        dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0] + m[3][0];
        dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1] + m[3][1];
        dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2] + m[3][2];
        w = src.x * m[0][3] + src.y * m[1][3] + src.z * m[2][3] + m[3][3];
        if (w != 1 && w != 0)
        {
            dst.x = dst.x / w;
            dst.y = dst.y / w;
            dst.z = dst.z / w;
        }
    }

    //transform a vector using the matrix
    template<typename S>
    void multDirMatrix(const Vec3<S> &src, Vec3<S> &dst) const
    {
        dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0];
        dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1];
        dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2];
    }

    //copy and transpose a the matrix
    Matrix44 transposed() const
    {
        Matrix44 transpMat;
        for( uint8_t i = 0; i < 4; ++i )
        {
            for( uint8_t j = 0; j < 4; ++j )
            {
                transpMat[i][j] = m[j][i];
            }
        }

        return transpMat;
    }

    //transpose the matrix
    Matrix44& transpose()
    {
        Matrix44 tmp(m[0][0],
                     m[1][0],
                     m[2][0],
                     m[3][0],
                     m[0][1],
                     m[1][1],
                     m[2][1],
                     m[3][1],
                     m[0][2],
                     m[1][2],
                     m[2][2],
                     m[3][2],
                     m[0][3],
                     m[1][3],
                     m[2][3],
                     m[3][3]);
        *this = tmp;

        return *this;
    }

    //get the inverse of the matrix (if there is one)
    Matrix44 inverse()
    {
        int i, j, k;
        Matrix44 s;
        Matrix44 t (*this);

        // Forward elimination
        for (i = 0; i < 3 ; i++)
        {
            int pivot = i;

            T pivotsize = t[i][i];

            if (pivotsize < 0)
                pivotsize = -pivotsize;

                for (j = i + 1; j < 4; j++) {
                    T tmp = t[j][i];

                    if (tmp < 0)
                        tmp = -tmp;

                        if (tmp > pivotsize) {
                            pivot = j;
                            pivotsize = tmp;
                        }
                }

            if (pivotsize == 0)
            {
                // Cannot invert singular matrix
                return Matrix44();
            }

            if (pivot != i)
            {
                for (j = 0; j < 4; j++)
                {
                    T tmp;

                    tmp = t[i][j];
                    t[i][j] = t[pivot][j];
                    t[pivot][j] = tmp;

                    tmp = s[i][j];
                    s[i][j] = s[pivot][j];
                    s[pivot][j] = tmp;
                }
            }

            for (j = i + 1; j < 4; j++)
            {
                T f = t[j][i] / t[i][i];

                for (k = 0; k < 4; k++)
                {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }

        // Backward substitution
        for (i = 3; i >= 0; --i)
        {
            T f;

            if ((f = t[i][i]) == 0)
            {
                // Cannot invert singular matrix
                return Matrix44();
            }

            for (j = 0; j < 4; j++)
            {
                t[i][j] /= f;
                s[i][j] /= f;
            }

            for (j = 0; j < i; j++)
            {
                f = t[j][i];

                for (k = 0; k < 4; k++) {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }

        return s;
    }

    //invert the matrix
    const Matrix44<T>& invert()
    {
        *this = inverse();
        return *this;
    }

    //output matrix
    friend std::ostream& operator << (std::ostream &s, const Matrix44 &x)
    {
        std::ios_base::fmtflags oldFlags = s.flags();
        int width = 12; // total with of the displayed number
        s.precision(5); // control the number of displayed decimals
        s.setf (std::ios_base::fixed);

        s << "(" << std::setw (width) << x[0][0] <<
             " " << std::setw (width) << x[0][1] <<
             " " << std::setw (width) << x[0][2] <<
             " " << std::setw (width) << x[0][3] << "\n" <<

             " " << std::setw (width) << x[1][0] <<
             " " << std::setw (width) << x[1][1] <<
             " " << std::setw (width) << x[1][2] <<
             " " << std::setw (width) << x[1][3] << "\n" <<

             " " << std::setw (width) << x[2][0] <<
             " " << std::setw (width) << x[2][1] <<
             " " << std::setw (width) << x[2][2] <<
             " " << std::setw (width) << x[2][3] << "\n" <<

             " " << std::setw (width) << x[3][0] <<
             " " << std::setw (width) << x[3][1] <<
             " " << std::setw (width) << x[3][2] <<
             " " << std::setw (width) << x[3][3] << ")\n";

        s.flags (oldFlags);
        return s;
    }
};

//shortcut for floating point matrix
typedef Matrix44<float> Matrix44f;
