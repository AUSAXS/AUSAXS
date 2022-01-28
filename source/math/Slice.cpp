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

std::string Slice::to_string() const {return operator Vector().to_string();}

//*** CONST & MUTABLE SLICE ***//
MutableSlice::MutableSlice(vector<double>& data, int N, int M, int start, int step, int length) 
    : Slice(N, M, start, step, length), data(data) {}
ConstSlice::ConstSlice(const vector<double>& data, int N, int M, int start, int step, int length) 
    : Slice(N, M, start, step, length), data(data) {}
VectorSlice::VectorSlice(const vector<double>& data) : Slice(data.size(), 1, 0, 1, data.size()), data(data) {}

double& MutableSlice::operator[](int i) {return data[start + i*step];}
const double& MutableSlice::operator[](int i) const {return data[start + i*step];}
const double& ConstSlice::operator[](int i) const {return data[start + i*step];}

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
Row::Row(vector<double>& data, int N, int M, int row) : MutableSlice(data, N, 1, row*M, 1, M) {}
ConstRow::ConstRow(const vector<double>& data, int N, int M, int row) : ConstSlice(data, N, 1, row*M, 1, M) {}
Column::Column(vector<double>& data, int N, int M, int col) : MutableSlice(data, 1, M, col, M, N) {}
ConstColumn::ConstColumn(const vector<double>& data, int N, int M, int col) : ConstSlice(data, 1, M, col, M, N) {}

Vector operator+(const Slice& s, const Vector& v) {return s.operator Vector().operator+(v);}
Vector operator-(const Slice& s, const Vector& v) {return s.operator Vector().operator-(v);}
Vector operator*(const Slice& s, const double a) {return s.operator Vector().operator*(a);}
Vector operator/(const Slice& s, const double a) {return s.operator Vector().operator/(a);}
bool operator==(const Slice& s, const Vector& v) {return s.operator Vector().operator==(v);}
bool operator!=(const Slice& s, const Vector& v) {return s.operator Vector().operator!=(v);}