#include "math/VectorN.h"
#include <algorithm>
#include <vector>

VectorN VectorN::operator+(const VectorN& w) const {
    std::vector<double> v(N);
    for (int i = 0; i < N; i++) {v[i] = data[i] + w[i];}
    return v;
}

VectorN& VectorN::operator+=(const VectorN& w) {
    for (int i = 0; i < N; i++) {data[i] += w[i];}
    return *this;
}

VectorN VectorN::operator-(const VectorN& w) const {
    std::vector<double> v(N); 
    std::transform(begin(), end(), w.begin(), std::back_inserter(v), [](double x1, double x2) {return x1-x2;});
    return v;
}

VectorN VectorN::operator*(const double& a) const {
    std::vector<double> v(N);
    for (int i = 0; i < N; i++) {v[i] = data[i]*a;}
    return v;
}

VectorN VectorN::operator/(const double& a) const {
    std::vector<double> v(N);
    for (int i = 0; i < N; i++) {v[i] = data[i]/a;}
    return v;
}