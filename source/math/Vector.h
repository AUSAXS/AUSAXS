#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>
#include <stdexcept>

// A basic vector class. Sizes are checked before each operation, so an std::invalid_argument is thrown if they do not match.
class Vector {
public:
    Vector(std::initializer_list<double> l) : data(l), N(l.size()) {}
    Vector(std::vector<double> w) : data(w), N(w.size()) {}
    Vector(const int& N) : N(N), data(N) {}
    ~Vector() {}

    // Assignment operator, w = v
    Vector& operator=(const Vector& v) {
        compatibility_check(v);
        data.assign(v.begin(), v.end());
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

    Vector operator-() const {
        Vector w(N);
        std::transform(begin(), end(), w.begin(), std::negate<double>());
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

    // Read-only indexing, w[i]
    const double& operator[](const int& i) const {return data[i];}
    
    // Read/write indexing, w[i] = ...
    double& operator[](const int& i) {return data[i];}

    // Approximate equality, w ~ v
    bool operator==(const Vector& v) const {
        compatibility_check(v);
        Vector a = operator-(v); // difference vector
        return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
    }

    // Approximate inequality operator, w != v
    bool operator!=(const Vector& v) const {return !operator==(v);}

    // Dot product
    double dot(const Vector& v) const {return std::inner_product(begin(), end(), v.begin(), 0);}

    // Norm (magnitude)
    double norm() const {return dot(*this);}

    // Euclidian distance to other vector
    double distance(const Vector& v) const {return sqrt(distance2(v));};

    // Euclidian distance squared to other vector
    double distance2(const Vector& v) const {
        compatibility_check(v);
        Vector w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), [] (double x1, double x2) {return pow((x1-x2), 2);});
        return std::accumulate(w.begin(), w.end(), 0);
    }

    // read-only iterators
    const std::vector<double>::const_iterator begin() const {return data.begin();}
    const std::vector<double>::const_iterator end() const {return data.end();}

    // read-write iterators
    std::vector<double>::iterator begin() {return data.begin();}
    std::vector<double>::iterator end() {return data.end();}

    const int size() const {return N;};

    // check if the vector is compatible with ours
    virtual void compatibility_check(const Vector& v) const {
        if (__builtin_expect(N != v.N, false)) {throw std::invalid_argument("Vector dimensions do not match.");}
    }

    const int N;

protected:
    std::vector<double> data;
    static constexpr double precision = 1e-9;
};