#include "Slice.h"
#include "Matrix.h"
#include "Vector.h"

//*** SLICE ***//
double Slice::dot(const Slice& s) {return dot(s.operator Vector());}
double Slice::dot(const Vector& v) {return Vector().dot(v);}
double Slice::norm() {return Vector().norm();}
Slice::operator Vector() const {
    Vector v(N);
    for (int i = 0; i < length; i++) {
        v[i] = operator[](i);
    }
    return v;
}

//*** CONST & MUTABLE SLICE ***//
MutableSlice::MutableSlice(Matrix* const parent, const int& start, const int& step, const int& length) : Slice(parent->N, parent->M, start, step, length), parent(parent) {}
ConstSlice::ConstSlice(const Matrix* parent, const int& start, const int& step, const int& length) : Slice(parent->N, parent->M, start, step, length), parent(parent) {}
VectorSlice::VectorSlice(const Vector* parent) : Slice(parent->N, 1, 0, 1, parent->N), parent(parent) {}

double& MutableSlice::operator[](const int& j) {return parent->_data[start + j*step];}
const double& MutableSlice::operator[](const int& j) const {return parent->data[start + j*step];}
const double& ConstSlice::operator[](const int& j) const {return parent->data[start + j*step];}

MutableSlice& MutableSlice::operator=(const Vector& v) {
    for (int i = 0; i < length; i++) {
        operator[](i) = v[i];
    }
    return *this;
}
MutableSlice& MutableSlice::operator+=(const Slice& s) {return operator+=(s.operator Vector());}
MutableSlice& MutableSlice::operator+=(const Vector& v) {
    for (int i = 0; i < length; i++) {
        operator[](i) += v[i];
    }
    return *this;
}

MutableSlice& MutableSlice::operator-=(const Slice& s) {return operator-=(s.operator Vector());}
MutableSlice& MutableSlice::operator-=(const Vector& v) {
    for (int i = 0; i < length; i++) {
        operator[](i) += v[i];
    }
    return *this;
}

//*** DERIVATIVES ***//
Row::Row(Matrix* const parent, const int& row) : MutableSlice(parent, row*parent->N, 1, parent->N) {}
Column::Column(Matrix* const parent, const int& col) : MutableSlice(parent, col, parent->M, parent->M) {}
ConstRow::ConstRow(const Matrix* parent, const int& row) : ConstSlice(parent, row*parent->N, 1, parent->N) {}
ConstColumn::ConstColumn(Matrix* const parent, const int& col) : ConstSlice(parent, col, parent->M, parent->M) {}

Vector operator+(const Slice& s, const Vector& v) {return s.operator Vector().operator+(v);}
Vector operator-(const Slice& s, const Vector& v) {return s.operator Vector().operator-(v);}
Vector operator*(const Slice& s, const double& a) {return s.operator Vector().operator*(a);}
Vector operator/(const Slice& s, const double& a) {return s.operator Vector().operator/(a);}
bool operator==(const Slice& s, const Vector& v) {return s.operator Vector().operator==(v);}
bool operator!=(const Slice& s, const Vector& v) {return s.operator Vector().operator!=(v);}