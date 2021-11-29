#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <stdexcept>

class Matrix {
public: 
    Matrix(std::initializer_list<std::initializer_list<double>> l) : N(N), M(M) {
        data.assign(l.begin(), l.end());
        for (const auto& e : data) {
            if (__builtin_expect(e.size() != data[0].size(), false)) {throw std::invalid_argument("Matrix columns must be of equal size!");}
        }
    }
    Matrix(const int& N, const int& M) : N(N), M(M), data(N, std::vector<double>(M, 0)) {}
    ~Matrix() {}

    // Assignment operator, B = A
    Matrix& operator=(const Matrix& A) {
        compatibility_check(A);
        data.assign(A.begin(), A.end());
        return *this;
    }

    // Plus operator, w + v
    Matrix operator+(const Matrix& A) const {
        compatibility_check(A);
        Matrix w(N, M);
        std::transform(begin(), end(), A.begin(), w.begin(), std::plus<double>());
        return w;
    }

    // Minus operator, w - v
    // Vector operator-(const Vector& v) const {
    //     compatibility_check(v);
    //     Vector w(N);
    //     std::transform(begin(), end(), v.begin(), w.begin(), std::minus<double>());
    //     return w;
    // }

    // Vector operator-() const {
    //     Vector w(N);
    //     std::transform(begin(), end(), w.begin(), std::negate<double>());
    //     return w;
    // }

    // // Scalar multiplication, w*a
    // Vector operator*(const double& a) const {
    //     Vector w(N);
    //     std::transform(begin(), end(), w.begin(), [&a](double x) {return x*a;});
    //     return w;
    // }

    // // Scalar division, w/a
    // Vector operator/(const double& a) const {
    //     Vector w(N);
    //     std::transform(begin(), end(), w.begin(), [&a](double x) {return x/a;});
    //     return w;
    // }

    // // Plus-assignment, w += v
    // Vector& operator+=(const Vector& v) {
    //     compatibility_check(v);
    //     std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
    //     return *this;
    // }

    // // Minus-assignment, w -= v
    // Vector& operator-=(const Vector& v) {
    //     compatibility_check(v);
    //     std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
    //     return *this;
    // }

    // // Read-only indexing, w[i]
    // const double& operator[](const int& i) const {return data[i];}
    
    // // Read/write indexing, w[i] = ...
    // double& operator[](const int& i) {return data[i];}

    // // Approximate equality, w ~ v
    // bool operator==(const Vector& v) const {
    //     compatibility_check(v);
    //     Vector a = operator-(v); // difference vector
    //     return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
    // }

    // // Approximate inequality operator, w != v
    // bool operator!=(const Vector& v) const {return !operator==(v);}

    // check if the matrix is compatible with ours
    virtual void compatibility_check(const Matrix& A) const {
        if (__builtin_expect(N != A.N || M != A.M, false)) {throw std::invalid_argument("Vector dimensions do not match.");}
    }

    // read-only iterators
    const std::vector<std::vector<double>>::const_iterator begin() const {return data.begin();}
    const std::vector<std::vector<double>>::const_iterator end() const {return data.end();}

    // read-write iterators
    std::vector<std::vector<double>>::iterator begin() {return data.begin();}
    std::vector<std::vector<double>>::iterator end() {return data.end();}

    const int N, M;

private: 
    std::vector<std::vector<double>> data;
    static constexpr double precision = 1e-9;
};