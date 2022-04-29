#pragma once

#include <initializer_list>

#include <math/Vector.h>

template<class T> class Matrix;

class Vector3 : public Vector<double> {
  public:
    /**
     * @brief Default constructor.
     */
    Vector3();

    /**
     * @brief Move constructor.
     */
    Vector3(Vector3&& v) noexcept;

    /**
     * @brief Copy constructor.
     */
    Vector3(const Vector3& v);

    /**
     * @brief List constructor.
     */
    Vector3(const std::initializer_list<double> l);

    /**
     * @brief Vector constructor.
     */
    Vector3(const Vector& v);

    /**
     * @brief Constructor. 
     * 
     * @param x The x-coordinate.
     * @param y The y-coordinate.
     * @param z The z-coordinate. 
     */
    Vector3(double x, double y, double z);

    /**
     * @brief Destructor.
     */
    ~Vector3() override;

    /**
     * @brief Set this vector equal to another. 
     */
    Vector3& operator=(const Vector3& v);

    /**
     * @brief Set this vector equal to an initializer list. 
     *        This allows the simple notation v = {1, 2, 3}.
     */
    Vector3& operator=(std::initializer_list<double> l);

    /**
     * @brief Get the Euclidian distance to another vector. 
     */
    double distance(const Vector3& v) const;

    /**
     * @brief Calculate the cross product of this vector with another. 
     */
    Vector3 cross(const Vector3& v) const;

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
    void rotate(const Matrix<double>& matrix);

    /**
     * @brief Output the string representation of this vector to a stream. 
     */
    friend std::ostream& operator<<(std::ostream& os, const Vector3& v) {os << v.to_string(); return os;}

    /**
     * @brief Check if this vector is approximately equal to another. 
     */
    bool operator==(const Vector3& v) const;

    /**
     * @brief Check if this vector is compatible with another vector. Two 3D vectors are always compatible. 
     */
    void compatibility_check(const Vector3&) const;

    /**
     * @brief Normalize this vector to unit length.
     */
    Vector3 normalize();

    /**
     * @brief Get a copy of the normalized form of this vector.
     */
    Vector3 normalize() const;

    /**
     * @brief Get a copy of the normalized form of this vector.
     */
    Vector3 normalize_copy() const;

    /**
     * @brief Generate a complete 3D basis from a single basis vector. 
     *        Implementation based on the "frisvad" algorithm from https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf.
     * 
     * @param n The first basis vector. 
     */
    static std::tuple<Vector3, Vector3, Vector3> generate_basis(const Vector3& v);

    /**
     * @brief Generate a complete 3D basis from this vector. 
     *        Implementation based on the "frisvad" algorithm from https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf.
     */
    std::tuple<Vector3, Vector3, Vector3> generate_basis();

    inline double& x() {return data[0];}
    inline const double& x() const {return data[0];}

    inline double& y() {return data[1];}
    inline const double& y() const {return data[1];}

    inline double& z() {return data[2];}
    inline const double& z() const {return data[2];}
};