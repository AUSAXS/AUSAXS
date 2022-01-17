#pragma once

class Vector;

#include <vector>
#include <stdexcept>

using std::vector;

class Slice {
  public:
    Slice(const int& N, const int& M, const int& start, const int& step, const int& length) 
        : N(N), M(M), start(start), step(step), length(length) {}

    double dot(const Slice& s);
    double dot(const Vector& v);
    double norm();
    operator Vector() const;
    virtual const double& operator[](const int& j) const = 0;
    void print() const;

    const size_t N, M;
    const int start, step, length;
};

Vector operator+(const Slice& s, const Vector& v);
Vector operator-(const Slice& s, const Vector& v);
Vector operator*(const Slice& s, const double& a);
Vector operator/(const Slice& s, const double& a);
bool operator==(const Slice& s, const Vector& v);
bool operator!=(const Slice& s, const Vector& v);

class VectorSlice : public Slice {
  public: 
    VectorSlice(const vector<double>& data);
    virtual ~VectorSlice() {}

    double& operator[](const int& i);
    const vector<double>& data;
};

class MutableSlice : public Slice {
  public:
    MutableSlice(vector<double>& data, const int& N, const int& M, const int& start, const int& step, const int& length);
    virtual ~MutableSlice() {}

    double& operator[](const int& i);
    const double& operator[](const int& i) const;
    MutableSlice& operator=(const Vector& v);
    MutableSlice& operator-=(const Slice& v);
    MutableSlice& operator-=(const Vector& v);
    MutableSlice& operator+=(const Slice& v);
    MutableSlice& operator+=(const Vector& v);

    vector<double>& data;
};

class ConstSlice : public Slice {
  public:
    ConstSlice(const vector<double>& data, const int& N, const int& M, const int& start, const int& step, const int& length);
    virtual ~ConstSlice() {}

    const double& operator[](const int& i) const;

    const vector<double>& data;
};

class ConstRow : public ConstSlice {
  public: 
    ConstRow(const vector<double>& data, const int& N, const int& M, const int& row);
    ConstRow(const ConstSlice& s) : ConstSlice(std::move(s)) {}
};

class ConstColumn : public ConstSlice {
  public: 
    ConstColumn(const vector<double>& data, const int& N, const int& M, const int& col);
    ConstColumn(const ConstSlice& s) : ConstSlice(std::move(s)) {}
};

class Row : public MutableSlice {
  public: 
    Row(vector<double>& data, const int& N, const int& M, const int& row);
    Row(MutableSlice& s) : MutableSlice(std::move(s)) {}

    Row& operator=(const Vector& v) {MutableSlice s(*this); MutableSlice::operator=(v); return *this;}
    Row& operator=(const Row& s) {
        for (int i = 0; i < length; i++) {
            operator[](i) = s[i];
        }
        return *this;
    }

    Row& operator=(const ConstRow& s) {
        for (int i = 0; i < length; i++) {
            operator[](i) = s[i];
        }
        return *this;
    }
};

class Column : public MutableSlice {
  public: 
    Column(vector<double>& data, const int& N, const int& M, const int& col);
    Column(MutableSlice& s) : MutableSlice(std::move(s)) {}

    Column& operator=(const Vector& v) {MutableSlice s(*this); MutableSlice::operator=(v); return *this;}
    Column& operator=(const Column& s) {
        for (int i = 0; i < length; i++) {
            operator[](i) = s[i];
        }
        return *this;
    }

    Column& operator=(const ConstColumn& s) {
        for (int i = 0; i < length; i++) {
            operator[](i) = s[i];
        }
        return *this;
    }
};