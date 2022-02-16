#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "Slice.h"

// A basic vector class. Sizes are checked before each operation, so an std::invalid_argument is thrown if they do not match.
class Vector {
  public:
    /**
     * @brief Move constructor.
     */
    Vector(Vector&& v) noexcept;

    /**
     * @brief Copy constructor.
     */
    Vector(const Vector& v);

    /**
     * @brief Construct a vector based on an initializer list.
     */
    Vector(const std::initializer_list<double> l);

    /**
     * @brief Construct a vector based on a std::vector. 
     */
    Vector(const std::vector<double>& v);

    /**
     * @brief Construct an empty vector of a given size. 
     */
    Vector(int n);

    /**
     * @brief Default constructor.
     */
    Vector();

    /**
     * @brief Destructor. 
     */
    virtual ~Vector();

    // Assignment operator, w = v
    Vector& operator=(const Vector& v);

    // Slice assignment operator
    Vector& operator=(const Slice& s);

    // Initializer list assignment operator
    Vector& operator=(std::initializer_list<double> l);

    // Plus operator, w + v
    Vector operator+(const Vector& v) const;

    // Minus operator, w - v
    Vector operator-(const Vector& v) const;

    // Minus operator, -w
    Vector operator-() const;

    // Vector multiplication, w*v
    Vector operator*(const Vector& v) const;

    // Scalar multiplication, w*a
    Vector operator*(double a) const;

    friend Vector operator*(double a, const Vector& v) {return v*a;}

    // Scalar division, w/a
    Vector operator/(double a) const;

    // Plus-assignment, w += v
    Vector& operator+=(const Vector& v);

    // Minus-assignment, w -= v
    Vector& operator-=(const Vector& v);

    // Scalar division-assignment, w /= a
    Vector& operator/=(double a);

    // Scalar multiplication-assignment, w /= a
    Vector& operator*=(double a);

    // Conversion to std::vector
    operator vector<double>();

    // Read-only indexing, w[i]
    const double& operator[](int i) const;
    
    // Read/write indexing, w[i] = ...
    double& operator[](int i);

    // Approximate equality, w ~ v
    bool operator==(const Vector& v) const;

    // Approximate inequality operator, w != v
    bool operator!=(const Vector& v) const;

    /**
     * @brief Get the dot product with another Vector.
     */
    double dot(const Vector& v) const;

    /**
     * @brief Get the norm of this Vector.
     */
    double norm() const;

    /**
     * @brief Get the magnitude of this Vector.
     */
    double magnitude() const;

    /**
     * @brief Get the Euclidian distance to another Vector.
     */
    double distance(const Vector& v) const;

    /**
     * @brief Get the squared Euclidian distance to another Vector.
     */
    double distance2(const Vector& v) const;

    /**
     * @brief Get a copy of this Vector.
     */
    Vector copy() const;

    /**
     * @brief Get a string representation of this Vector.
     */
    std::string to_string(std::string message = "") const;

    /**
     * @brief Output the string representation of this Vector to a stream. 
     */
    friend std::ostream& operator<<(std::ostream& os, const Vector& v) {os << v.to_string(); return os;}

    // Read-only iterator
    const std::vector<double>::const_iterator begin() const;

    // Read-only iterator
    const std::vector<double>::const_iterator end() const;

    // Read-write iterator
    std::vector<double>::iterator begin();

    // Read-write iterator
    std::vector<double>::iterator end();

    /**
     * @brief Get the size of this Vector.
     */
    size_t size() const;

    /**
     * @brief Get the dimension of this Vector.
     */
    size_t dim() const;

    const size_t& N = _N; // read-only access to the dimension
    const std::vector<double>& data = _data; // read-only access to the data container

  protected:
    size_t _N;
    std::vector<double> _data;
    static constexpr double precision = 1e-9;

    // check if the vector is compatible with ours
    virtual void compatibility_check(const Vector& v) const;
};