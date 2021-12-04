#pragma once

#include "Slice.h"
#include "Matrix.h"
#include "Vector.h"

ConstRowSlice::ConstRowSlice(const Matrix* parent, const int& row) : Slice(parent), row(row) {}  
ConstColumnSlice::ConstColumnSlice(Matrix* const parent, const int& col) : Slice(parent), col(col) {}
ColumnSlice::ColumnSlice(Matrix* const parent, const int& col) : Slice(parent), col(col) {}
RowSlice::RowSlice(Matrix* const parent, const int& row) : Slice(parent), row(row) {}

RowSlice& RowSlice::operator=(const RowSlice& s) {operator=(Vector(s));}
RowSlice& RowSlice::operator=(const Vector& v) {
    if (__builtin_expect(N != v.N, false)) {throw std::invalid_argument("Vector dimension do not match (got: [" + std::to_string(v.N) + ", expected: " + std::to_string(N) + ")");}
    Matrix p = *parent;
    for (int i = 0; i < M; i++) {
        p._data[row*M + i] = v[i];
    }
    return *this;
}

ColumnSlice& ColumnSlice::operator=(const ColumnSlice& s) {operator=(Vector(s));}
ColumnSlice& ColumnSlice::operator=(const Vector& v) {
    if (__builtin_expect(M != v.N, false)) {throw std::invalid_argument("Vector dimension do not match (got: [" + std::to_string(v.N) + ", expected: " + std::to_string(N) + ")");}
    Matrix p = *parent;
    for (int i = 0; i < N; i++) {
        p._data[i*M + col] = v[i];
    }
    return *this;
}

std::vector<double> RowSlice::copy() const {
    std::vector<double> v(N);
    Matrix p = *parent;
    for (int i = 0; i < N; i++) {
        v[i] = p.data[row*M + i];
    }
    return v;
}

std::vector<double> ConstRowSlice::copy() const {
    std::vector<double> v(N);
    Matrix p = *parent;
    for (int i = 0; i < N; i++) {
        v[i] = p.data[row*M + i];
    }
    return v;
}

std::vector<double> ColumnSlice::copy() const {
    std::vector<double> v(N);
    Matrix p = *parent;
    for (int i = 0; i < N; i++) {
        v[i] = p.data[i*M + col];
    }
    return v;
}

std::vector<double> ConstColumnSlice::copy() const {
    std::vector<double> v(N);
    Matrix p = *parent;
    for (int i = 0; i < N; i++) {
        v[i] = p.data[i*M + col];
    }
    return v;
}

RowSlice::operator Vector() const {return Vector(copy());}
ConstRowSlice::operator Vector() const {return Vector(copy());}
ColumnSlice::operator Vector() const {return Vector(copy());}
ConstColumnSlice::operator Vector() const {return Vector(copy());}

const double& ConstRowSlice::operator[](const int& j) const {return parent->data[row*M + j];}
const double& RowSlice::operator[](const int& j) const {return parent->data[row*M + j];}
double& RowSlice::operator[](const int& j) {return parent->_data[row*M + j];}

const double& ConstColumnSlice::operator[](const int& j) const {return parent->data[j*M + col];}
const double& ColumnSlice::operator[](const int& j) const {return parent->data[j*M + col];}
double& ColumnSlice::operator[](const int& j) {return parent->_data[j*M + col];}
