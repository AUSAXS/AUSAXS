#include <initializer_list>

template <class T>
class vector3 {
public:
    T x, y, z; 
    vector3(std::initializer_list<T> l) : x(l[0]), y(l[1]), z(l[2]) {}
    vector3(T x, T y, T z) : x(x), y(y), z(z) {}

    vector3 operator+(const vector3& w) const {return vector3(x+w.x, y+w.y, z+w.z);}
    vector3 operator-(const vector3& w) const {return vector3(x-w2.x, y-w.y, z-w.z);}
    vector3 operator*(const double& a) const {return vector3(a*x, a*y, a*z);}
    vector3 operator/(const double& a) const {return vector3(x/a, y/a, z/a);}
};