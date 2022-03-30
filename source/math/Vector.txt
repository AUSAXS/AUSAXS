#include <initializer_list>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "math/Slice.h"
#include "math/Vector.h"

Vector::Vector(Vector&& v) noexcept : _N(v._N), _data(std::move(v._data)) {}

Vector::Vector(const Vector& v) : _N(v.size()), _data(v._data) {}

Vector::Vector(const std::initializer_list<double> l) : _N(l.size()), _data(l) {}

Vector::Vector(const std::vector<double>& v) : _N(v.size()), _data(v) {}

Vector::Vector(int n) : _N(n), _data(n) {}

Vector::Vector() : _N(0), _data(0) {}

Vector::~Vector() = default;

Vector& Vector::operator=(const Vector& v) {
    _N = v.N;
    _data = v.data;
    return *this;
}

Vector& Vector::operator=(const Slice& s) {
    if (__builtin_expect(!(s.N == 1 || s.M == 1), false)) {throw std::invalid_argument("Only 1D slices can be assigned to vectors. Size: " + std::to_string(s.N) + ", " + std::to_string(s.M));}
    _N = std::max(s.N, s.M);
    _data = std::vector<double>(N);
    for (size_t i = 0; i < N; i++) {
        _data[i] = s[i];
    }
    return *this;
}

Vector& Vector::operator=(std::initializer_list<double> l) {
    _N = l.size();
    _data.assign(l);
    return *this;
}

Vector Vector::operator+(const Vector& v) const {
    compatibility_check(v);
    Vector w(N);
    std::transform(begin(), end(), v.begin(), w.begin(), std::plus<double>());
    return w;
}

Vector Vector::operator-(const Vector& v) const {
    compatibility_check(v);
    Vector w(N);
    std::transform(begin(), end(), v.begin(), w.begin(), std::minus<double>());
    return w;
}

Vector Vector::operator-() const {
    Vector w(N);
    std::transform(begin(), end(), w.begin(), std::negate<double>());
    return w;
}

Vector Vector::operator*(const Vector& v) const {
    compatibility_check(v);
    Vector w(N);
    std::transform(begin(), end(), v.begin(), w.begin(), std::multiplies<double>());
    return w;
}

Vector Vector::operator*(double a) const {
    Vector w(N);
    std::transform(begin(), end(), w.begin(), [&a](double x) {return x*a;});
    return w;
}

Vector Vector::operator/(double a) const {
    Vector w(N);
    std::transform(begin(), end(), w.begin(), [&a](double x) {return x/a;});
    return w;
}

Vector& Vector::operator+=(const Vector& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
    return *this;
}

Vector& Vector::operator-=(const Vector& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
    return *this;
}

Vector& Vector::operator/=(double a) {
    std::transform(begin(), end(), begin(), [&a](double x) {return x/a;});
    return *this;
}

Vector& Vector::operator*=(double a) {
    std::transform(begin(), end(), begin(), [&a](double x) {return x*a;});
    return *this;
}

Vector::operator vector<double>() {
    return data;
}

const double& Vector::operator[](int i) const {return data[i];}

double& Vector::operator[](int i) {return _data[i];}

bool Vector::operator==(const Vector& v) const {
    compatibility_check(v);
    Vector a = operator-(v); // difference vector
    return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
}

bool Vector::operator!=(const Vector& v) const {return !operator==(v);}

double Vector::dot(const Vector& v) const {
    compatibility_check(v);
    return std::inner_product(begin(), end(), v.begin(), 0.0);
}

double Vector::norm() const {return sqrt(dot(*this));}

double Vector::magnitude() const {return norm();}

double Vector::distance(const Vector& v) const {return sqrt(distance2(v));};

double Vector::distance2(const Vector& v) const {
    compatibility_check(v);
    Vector w(N);
    std::transform(begin(), end(), v.begin(), w.begin(), [] (double x1, double x2) {return pow((x1-x2), 2);});
    return std::accumulate(w.begin(), w.end(), 0);
}

Vector Vector::copy() const {
    Vector w(N);
    std::copy(begin(), end(), w.begin());
    return w;
}

std::string Vector::to_string(std::string message) const {
    if (message != "") {std::cout << message << std::endl;}
    std::stringstream s; s << "( ";
    for (const auto& e : data) {
        s << std::setprecision(8) << e << " ";
    }
    s << ")";
    return s.str();
}

const std::vector<double>::const_iterator Vector::begin() const {return data.begin();}

const std::vector<double>::const_iterator Vector::end() const {return data.end();}

std::vector<double>::iterator Vector::begin() {return _data.begin();}

std::vector<double>::iterator Vector::end() {return _data.end();}

size_t Vector::size() const {return N;};

size_t Vector::dim() const {return size();}

void Vector::compatibility_check(const Vector& v) const {
    if (__builtin_expect(N != v.N, false)) {
        throw std::invalid_argument("Vector dimensions do not match (got: " + std::to_string(v.N) + ", expected: " + std::to_string(N) + ").");
    }
}