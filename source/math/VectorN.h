#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>

class VectorN {
public:
    VectorN(std::initializer_list<double> l) : data(l), N(l.size()) {}
    VectorN(std::vector<double> w) : data(w), N(w.size()) {}
    VectorN(const int& N) : N(N) {data = std::vector<double>(N);}

    VectorN operator+(const VectorN& v) const {
        std::vector<double> w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), std::plus<double>());
        return VectorN(w);
    }

    VectorN operator-(const VectorN& v) const {
        std::vector<double> w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), std::minus<double>());
        return VectorN(w);
    }

    VectorN operator*(const double& a) const {
        std::vector<double> w(N);
        std::transform(begin(), end(), w.begin(), [&a](double x) {return x*a;});
        return VectorN(w);
    }

    VectorN operator/(const double& a) const {
        std::vector<double> w(N);
        std::transform(begin(), end(), w.begin(), [&a](double x) {return x/a;});
        return VectorN(w);
    }

    VectorN& operator+=(const VectorN& v) {
        std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
        return *this;
    }

    VectorN& operator-=(const VectorN& v) {
        std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
        return *this;
    }

    const double& operator[](const int& i) const {return data[i];}
    double& operator[](const int& i) {return data[i];}

    bool operator==(const VectorN& v) { // returns true if the vectors are APPROXIMATELY equal
        VectorN a = operator-(v); // difference vector
        return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
    }

    double dot(const VectorN& v) const {return std::inner_product(begin(), end(), v.begin(), 0);}
    double norm() const {return dot(*this);}
    double distance(const VectorN& v) const {return sqrt(distance2(v));};
    double distance2(const VectorN& v) const {
        std::vector<double>	w(N);
        std::transform(begin(), end(), v.begin(), w.begin(), [] (double x1, double x2) {return pow((x1-x2), 2);});
        return std::accumulate(w.begin(), w.end(), 0);
    }

    const std::vector<double>::const_iterator begin() const {return data.begin();}
    const std::vector<double>::const_iterator end() const {return data.end();}
    std::vector<double>::iterator begin() {return data.begin();}
    std::vector<double>::iterator end() {return data.end();}

    const int size() const {return N;};

protected:
    const int N;
    std::vector<double> data;
    static constexpr double precision = 1e-9;
};