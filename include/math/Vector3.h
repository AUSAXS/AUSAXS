#pragma once

#include <initializer_list>
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
     * @brief Calculate the Euclidian distance to another vector. 
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
    void rotate(Vector3& axis, const double& angle);

    /**
     * @brief Rotate this vector by a rotation matrix.
     * 
     * @param matrix The rotation matrix. 
     */
    void rotate(Matrix& matrix);

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
     * @brief Get a string representation of this vector. 
     */
    std::string to_string() const {return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";}

    /**
     * @brief Normalize this vector to unit length.
     */
    void normalize() {
        operator=(operator/(norm()));
    }

    /**
     * @brief Normalize this vector to unit length.
     */
    Vector3 normalize() const {
        return operator/(norm());
    }

    // Allow mutable access to the data through the simple v.x, v.y, and v.z notation. 
    double& x = _data[0];
    double& y = _data[1];
    double& z = _data[2];
};