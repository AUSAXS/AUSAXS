#include <initializer_list>
#include "math/VectorN.h"

template <class T>
class Vector3 : public VectorN {
public:
    T x, y, z; 
    Vector3(std::initializer_list<T> l) : x(l[0]), y(l[1]), z(l[2]) {}
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

    Vector3 operator+(const Vector3& w) const {return vector3(x+w.x, y+w.y, z+w.z);}
    Vector3 &operator+=(const Vector3& w) {x+=w.x, y+=w.y, z+=w.z; return *this;}
    Vector3 operator-(const Vector3& w) const {return vector3(x-w2.x, y-w.y, z-w.z);}
    Vector3 &operator-=(const Vector3& w) {x-=w.x, y-=w.y, z-=w.z; return *this;}
    Vector3 operator*(const double& a) const {return vector3(a*x, a*y, a*z);}
    Vector3 operator/(const double& a) const {return vector3(x/a, y/a, z/a);}
};