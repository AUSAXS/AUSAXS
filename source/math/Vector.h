#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "Slice.h"

// A basic vector class. Sizes are checked before each operation, so an std::invalid_argument is thrown if they do not match.
class Vector {
public:
    Vector(const Vector& v) : _N(v.size()), _data(v._data) {} // copy constructor
    Vector(const std::initializer_list<double> l) : _N(l.size()), _data(l) {} // initializer list {a, b, c, d}
    Vector(const std::vector<double>& v) : _N(v.size()), _data(v) {} // std::vector --> Vector constructor
    Vector(const int& n) : _N(n), _data(n) {} // dimensional constructor
    Vector() : _N(0), _data(0) {} // default constructor
    virtual ~Vector() {}

    // Assignment operator, w = v
    Vector& operator=(const Vector& v) {
        _N = v.N;
        _data = v.data;
        return *this;
    }

    template<class T>
    Vector& operator=(const Slice& s) {
        if (__builtin_expect(!(s.N == 1 || s.M == 1), false)) {throw std::invalid_argument("Slice is not convertible to a vector.");}
        _N = std::max(s.N, s.M);
        _data = std::vector<double>(N);
        for (size_t i = 0; i < N; i++) {
            _data[i] = s[i];
        }
        return *this;
    }

    Vector& operator=(std::initializer_list<double> l) {
        _N = l.size();
        _data.assign(l);
        return *this;
    }

    // Plus operator, w + v
    Vector operator+(const Vector& v) const {
        compatibility_check(v);
        Vector w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), std::plus<double>());
        return w;
    }

    // Minus operator, w - v
    Vector operator-(const Vector& v) const {
        compatibility_check(v);
        Vector w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), std::minus<double>());
        return w;
    }

    // Minus operator, -w
    Vector operator-() const {
        Vector w(N);
        std::transform(begin(), end(), w.begin(), std::negate<double>());
        return w;
    }

    // Vector multiplication, w*v
    Vector operator*(const Vector& v) const {
        compatibility_check(v);
        Vector w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), std::multiplies<double>());
        return w;
    }

    // Scalar multiplication, w*a
    Vector operator*(const double& a) const {
        Vector w(N);
        std::transform(begin(), end(), w.begin(), [&a](double x) {return x*a;});
        return w;
    }

    // Scalar division, w/a
    Vector operator/(const double& a) const {
        Vector w(N);
        std::transform(begin(), end(), w.begin(), [&a](double x) {return x/a;});
        return w;
    }

    // Plus-assignment, w += v
    Vector& operator+=(const Vector& v) {
        compatibility_check(v);
        std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
        return *this;
    }

    // Minus-assignment, w -= v
    Vector& operator-=(const Vector& v) {
        compatibility_check(v);
        std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
        return *this;
    }

    // Scalar division-assignment, w /= a
    Vector& operator/=(const double& a) {
        std::transform(begin(), end(), begin(), [&a](double x) {return x/a;});
        return *this;
    }

    // Scalar multiplication-assignment, w /= a
    Vector& operator*=(const double& a) {
        std::transform(begin(), end(), begin(), [&a](double x) {return x*a;});
        return *this;
    }

    // Read-only indexing, w[i]
    const double& operator[](const int& i) const {return data[i];}
    
    // Read/write indexing, w[i] = ...
    double& operator[](const int& i) {return _data[i];}

    // Approximate equality, w ~ v
    bool operator==(const Vector& v) const {
        compatibility_check(v);
        Vector a = operator-(v); // difference vector
        return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
    }

    // Approximate inequality operator, w != v
    bool operator!=(const Vector& v) const {return !operator==(v);}

    // Dot product
    double dot(const Vector& v) const {
        compatibility_check(v);
        return std::inner_product(begin(), end(), v.begin(), 0.0);
    }

    // Norm (magnitude)
    double norm() const {return sqrt(dot(*this));}

    // Euclidian distance to other vector
    double distance(const Vector& v) const {return sqrt(distance2(v));};

    // Euclidian distance squared to other vector
    double distance2(const Vector& v) const {
        compatibility_check(v);
        Vector w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), [] (double x1, double x2) {return pow((x1-x2), 2);});
        return std::accumulate(w.begin(), w.end(), 0);
    }

    // Returns a copy of this vector
    Vector copy() const {
        Vector w(N);
        std::copy(begin(), end(), w.begin());
        return w;
    }

    // Print this vector to the terminal
    void print(const std::string& message = "") const {
        if (message != "") {std::cout << message << std::endl;}
        for (const auto& e : data) {
            std::cout << std::setw(8) << e << " ";
        }
        std::cout << std::endl;
    }

    // read-only iterators
    const std::vector<double>::const_iterator begin() const {return data.begin();}
    const std::vector<double>::const_iterator end() const {return data.end();}

    // read-write iterators
    std::vector<double>::iterator begin() {return _data.begin();}
    std::vector<double>::iterator end() {return _data.end();}

    size_t size() const {return N;};

    const size_t& N = _N; // read-only access to the dimension
    const std::vector<double>& data = _data; // read-only access to the data container

protected:
    size_t _N;
    std::vector<double> _data;
    static constexpr double precision = 1e-9;

    // check if the vector is compatible with ours
    virtual void compatibility_check(const Vector& v) const {
        if (__builtin_expect(N != v.N, false)) {
            throw std::invalid_argument("Vector dimensions do not match (got: " + std::to_string(v.N) + ", expected: " + std::to_string(N) + ").");
        }
    }
};