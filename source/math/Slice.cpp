#include "math/Slice.h"
#include "math/Matrix.h"
#include "math/Vector.h"

#include <vector>

using std::vector;

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
MutableSlice::MutableSlice(vector<double>& data, const int& N, const int& M, const int& start, const int& step, const int& length) 
    : Slice(N, M, start, step, length), data(data) {}
ConstSlice::ConstSlice(const vector<double>& data, const int& N, const int& M, const int& start, const int& step, const int& length) 
    : Slice(N, M, start, step, length), data(data) {}
VectorSlice::VectorSlice(const vector<double>& data) : Slice(data.size(), 1, 0, 1, data.size()), data(data) {}

double& MutableSlice::operator[](const int& i) {return data[start + i*step];}
const double& MutableSlice::operator[](const int& i) const {return data[start + i*step];}
const double& ConstSlice::operator[](const int& i) const {return data[start + i*step];}

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
Row::Row(vector<double>& data, const int& N, const int& M, const int& row) : MutableSlice(data, N, M, row*M, 1, M) {}
ConstRow::ConstRow(const vector<double>& data, const int& N, const int& M, const int& row) : ConstSlice(data, N, M, row*M, 1, M) {}
Column::Column(vector<double>& data, const int& N, const int& M, const int& col) : MutableSlice(data, N, M, col, M, N) {}
ConstColumn::ConstColumn(const vector<double>& data, const int& N, const int& M, const int& col) : ConstSlice(data, N, M, col, M, N) {}

Vector operator+(const Slice& s, const Vector& v) {return s.operator Vector().operator+(v);}
Vector operator-(const Slice& s, const Vector& v) {return s.operator Vector().operator-(v);}
Vector operator*(const Slice& s, const double& a) {return s.operator Vector().operator*(a);}
Vector operator/(const Slice& s, const double& a) {return s.operator Vector().operator/(a);}
bool operator==(const Slice& s, const Vector& v) {return s.operator Vector().operator==(v);}
bool operator!=(const Slice& s, const Vector& v) {return s.operator Vector().operator!=(v);}