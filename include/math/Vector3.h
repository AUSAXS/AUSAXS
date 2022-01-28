#pragma once

#include <initializer_list>
#include <tuple>

#include "math/Vector.h"
#include "math/Matrix.h"

class Vector3 : public Vector {
  public:
    Vector3(Vector3&& v) noexcept : Vector(std::move(v)) {} // move constructor
    Vector3(const Vector3& v) : Vector(v.data) {} // copy constructor
    Vector3() : Vector(3) {} // default empty constructor
    Vector3(const std::initializer_list<double> l) : Vector(l) {} // initializer list {a, b, c, d}
    Vector3(const Vector& v) : Vector(v) {}
    Vector3(double x, double y, double z) : Vector({x, y, z}) {}
    ~Vector3() override {}

    /**
     * @brief Set this vector equal to another. 
     */
    Vector3& operator=(const Vector3& v) {
        _data = v.data;
        return *this;
    }

    /**
     * @brief Set this vector equal to an initializer list. 
     *        This allows the simple notation v = {1, 2, 3}.
     */
    Vector3& operator=(std::initializer_list<double> l) {
        _N = l.size();
        _data.assign(l);
        return *this;
    }

    /**
     * @brief Get the Euclidian distance to another vector. 
     */
    double distance(const Vector3& v) const {return sqrt(pow(x-v.x, 2) + pow(y-v.y, 2) + pow(z-v.z, 2));}

    /**
     * @brief Calculate the cross product of this vector with another. 
     */
    Vector3 cross(const Vector3& v) const {return {y*v.z - v.y*z, z*v.x - v.z*x, x*v.y - v.x*y};}

    /**
     * @brief Rotate this vector around an axis. 
     * 
     * @param axis The rotation axis. 
     * @param angle The angle to rotate. 
     */
    void rotate(const Vector3& axis, double angle);

    /**
     * @brief Rotate this vector by a rotation matrix.
     * 
     * @param matrix The rotation matrix. 
     */
    void rotate(const Matrix& matrix);

    // /**
    //  * @brief Get a string representation of this vector. 
    //  */
    // std::string to_string() const {return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";}

    /**
     * @brief Output the string representation of this vector to a stream. 
     */
    friend std::ostream& operator<<(std::ostream& os, const Vector3& v) {os << v.to_string(); return os;}

    /**
     * @brief Check if this vector is approximately equal to another. 
     */
    bool operator==(const Vector3& v) const {return abs(x-v.x) + abs(y-v.y) + abs(z-v.z) < precision; }

    /**
     * @brief Check if this vector is compatible with another vector. Two 3D vectors are always compatible. 
     */
    void compatibility_check(const Vector3&) const {}

    /**
     * @brief Normalize this vector to unit length.
     */
    Vector3 normalize() {
        operator=(operator/(norm()));
        return *this;
    }

    /**
     * @brief Get a copy of the normalized form of this vector.
     */
    Vector3 normalize() const {
        return operator/(norm());
    }

    /**
     * @brief Get a copy of the normalized form of this vector.
     */
    Vector3 normalize_copy() const {
        return operator/(norm());
    }

    /**
     * @brief Generate a complete 3D basis from a single basis vector. 
     *        Implementation based on the "frisvad" algorithm from https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf.
     * 
     * @param n The first basis vector. 
     */
    static std::tuple<Vector3, Vector3, Vector3> generate_basis(const Vector3& v) {
        Vector3 n = v.normalize_copy();

        // Handle the singularity
        if (n.z < -0.9999999) { 
            Vector3 b1(0, -1, 0);
            Vector3 b2(-1, 0, 0);
            return std::make_tuple(n, b1, b2);
        }
        const float a = 1/(1 + n.z);
        const float b = -n.x*n.y*a;
        Vector3 b1(1-n.x*n.x*a, b, -n.x);
        Vector3 b2(b, 1-n.y * n.y*a, -n.y);
        return std::make_tuple(n, b1, b2);
    }

    /**
     * @brief Generate a complete 3D basis from this vector. 
     *        Implementation based on the "frisvad" algorithm from https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf.
     */
    std::tuple<Vector3, Vector3, Vector3> generate_basis() {
        return Vector3::generate_basis(*this);
    }

    // Allow mutable access to the data through the simple v.x, v.y, and v.z notation. 
    double& x = _data[0];
    double& y = _data[1];
    double& z = _data[2];
};