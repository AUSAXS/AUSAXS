#include "Slice.h"
#include "Matrix.h"
#include "Vector.h"

//*** SLICE ***//
double Slice::dot(const Slice& s) {return dot(s.operator Vector());}
double Slice::dot(const Vector& v) {return operator Vector().dot(v);}
double Slice::norm() {return operator Vector().norm();}
Slice::operator Vector() const {
    Vector v(length);
    for (int i = 0; i < length; i++) {
        v[i] = operator[](i);
    }
    return v;
}

void Slice::print() const {return operator Vector().print();}

//*** CONST & MUTABLE SLICE ***//
MutableSlice::MutableSlice(Matrix* const parent, const int& start, const int& step, const int& length) : Slice(parent->N, parent->M, start, step, length), parent(parent) {}
ConstSlice::ConstSlice(const Matrix* parent, const int& start, const int& step, const int& length) : Slice(parent->N, parent->M, start, step, length), parent(parent) {}
VectorSlice::VectorSlice(const Vector* parent) : Slice(parent->N, 1, 0, 1, parent->N), parent(parent) {}

double& MutableSlice::operator[](const int& i) {return parent->_data[start + i*step];}
const double& MutableSlice::operator[](const int& i) const {return parent->data[start + i*step];}
const double& ConstSlice::operator[](const int& i) const {return parent->data[start + i*step];}

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
        operator[](i) -= v[i];
    }
    return *this;
}

//*** DERIVATIVES ***//
Row::Row(Matrix* const parent, const int& row) : MutableSlice(parent, row*parent->M, 1, parent->M) {}
ConstRow::ConstRow(const Matrix* parent, const int& row) : ConstSlice(parent, row*parent->M, 1, parent->M) {}
Column::Column(Matrix* const parent, const int& col) : MutableSlice(parent, col, parent->M, parent->N) {}
ConstColumn::ConstColumn(const Matrix* parent, const int& col) : ConstSlice(parent, col, parent->M, parent->N) {}

Vector operator+(const Slice& s, const Vector& v) {return s.operator Vector().operator+(v);}
Vector operator-(const Slice& s, const Vector& v) {return s.operator Vector().operator-(v);}
Vector operator*(const Slice& s, const double& a) {return s.operator Vector().operator*(a);}
Vector operator/(const Slice& s, const double& a) {return s.operator Vector().operator/(a);}
bool operator==(const Slice& s, const Vector& v) {return s.operator Vector().operator==(v);}
bool operator!=(const Slice& s, const Vector& v) {return s.operator Vector().operator!=(v);}