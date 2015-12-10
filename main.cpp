#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <random>
#include "geometry.h"

#define M_PI 3.14159265

using namespace std;

//infinity is the largest floating point number
const float kInfinity = numeric_limits<float>::max();

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0,1);

inline float clamp(const float &lo, const float &hi, const float &v)
{
    return max(lo, min(hi,v));
}

// convert degrees to radians
inline float deg2rad(const float &deg)
{
    return deg*M_PI/180;
}

inline Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{
    return a*(1 - mixValue) + b*mixValue;
}

// the options for the picture we will generate
struct Options
{
    uint32_t width;
    uint32_t height;
    float fov;
    Matrix44f cameraToWorld;
};

//the basic object class, from which other shapes are derived
class Object
{
public:
    //the empty constructor
    Object() : color(dis(gen), dis(gen), dis(gen)) {}

    //the destructor
    virtual ~Object()
    {

    }

    virtual bool intersect(const Vec3f &, const Vec3f &, float &) const = 0;
    virtual void getSurfaceData(const Vec3f &, Vec3f &, Vec2f &) const = 0;
    Vec3f color;
};

//finds the real roots of a quadratic equation
bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
    float discr = b*b - 4*a*c;
    if( discr < 0 )
    {
        return false;
    }
    else
    {
        x0 = (-1*b + sqrt(discr))/(2*a);
        x1 = (-1*b - sqrt(discr))/(2*a);
    }
    if( x0 > x1 )
    {
        swap(x0,x1);
    }

    return true;
}

class Sphere : public Object
{
public:
    //constructor
    Sphere(const Vec3f &c, const float &r)
    {
        radius = r;
        radius2 = r;
        center = c;
    }

    //specific intersection checker for spheres
    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const
    {
        float t0, t1;
        Vec3f L = center - orig;
        float tca = L.dotProduct(dir);
        if(tca < 0)
        {
            return false;
        }

        float d2 = L.dotProduct(L) - tca*tca;
        if( d2 > radius2)
        {
            return false;
        }

        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;

        if(t0 > t1)
        {
            swap(t0,t1);
        }

        if(t0 < 0)
        {
            t0 = t1;
            if( t0 < 0)
            {
                return false;
            }
        }

        t = t0;

        return true;
    }

    //find data about the surface of the sphere
    void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const
    {
        Nhit = Phit - center;
        Nhit.normalize();

        tex.x = (1 + atan2(Nhit.z, Nhit.x)/M_PI)*0.5;
        tex.y = acosf(Nhit.y)/M_PI;
    }

    float radius, radius2;
    Vec3f center;
};

//checks for ray intersecting object, gives nearest intersection point
bool trace(const Vec3f &orig, const Vec3f &dir, const vector<unique_ptr<Object>> &objects, float &tNear, const Object *&hitObject)
{
    tNear = kInfinity;
    vector<unique_ptr<Object>>::const_iterator iter = objects.begin();
    for(; iter != objects.end(); ++iter)
    {
        float t = kInfinity;
        if( (*iter)->intersect(orig, dir, t) && t < tNear)
        {
            hitObject = iter->get();
            tNear = t;
        }
    }

    return (hitObject != nullptr);
}

//casts a ray into the scene
Vec3f castRay(const Vec3f &orig, const Vec3f &dir, const vector<unique_ptr<Object>> &objects)
{
    Vec3f hitColor = 0;
    const Object *hitObject = nullptr;
    float t;
    if( trace(orig, dir, objects, t, hitObject) )
    {
        Vec3f Phit = orig + dir*t;
        Vec3f Nhit;
        Vec2f tex;
        hitObject->getSurfaceData(Phit, Nhit, tex);
        float scale = 4;
        float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y*scale, 1) > 0.5);
        hitColor = std::max(0.f, Nhit.dotProduct(-dir)) * mix(hitObject->color, 0.8*hitObject->color, pattern);
    }

    return hitColor;
}

void render(const Options &options, const vector<unique_ptr<Object>> &objects)
{
    Vec3f *framebuffer = new Vec3f[options.width * options.height];
    Vec3f *pix = framebuffer;
    float scale = tan(deg2rad(options.fov*0.5));
    float imageAspectRatio = options.width / (float)options.height;
    Vec3f orig;
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);

    //the double-nested for loop where each pixel is calculated, this is where parallelization will occur
    #pragma omp parallel
    {
        for(uint32_t j = 0; j < options.height; ++j)
        {
            #pragma omp for ordered schedule(dynamic)
            for(uint32_t i = 0; i <options.width; ++i)
            {
                float x = (2*(i + 0.5)/(float)options.width - 1)*imageAspectRatio*scale;
                float y = (1 - 2*(j + 0.5)/(float)options.height)*scale;

                Vec3f dir;
                options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
                dir.normalize();

                *(pix++) = castRay(orig, dir, objects);
            }
        }

        //save to a PPM file
        ofstream ofs("./out.ppm", ios::out | ios::binary);
        ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
        for (uint32_t i = 0; i < options.height*options.width; ++i)
        {
            char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
            char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
            char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
            ofs << r << g << b;
        }
        ofs.close();
    }
    delete [] framebuffer;
}
int main(int argc, char **argv)
{
    vector<unique_ptr<Object>> objects;

    //make random spheres
    uint32_t numSpheres = 128;
    gen.seed(0);
    #pragma omp for
    for(uint32_t i = 0; i < numSpheres; ++i)
    {
        Vec3f randPos((0.5 - dis(gen))*10, (0.5 - dis(gen))*10, (0.5 + dis(gen)*10));
        float randRadius = (0.5 + dis(gen)*0.5);
        objects.push_back(unique_ptr<Object>(new Sphere(randPos, randRadius)));
    }

    //set image options
    Options options;
    options.width = 1024;
    options.height = 768;
    options.fov = 60;
    options.cameraToWorld = Matrix44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

    //render the image
    render(options, objects);

    return 0;
}
