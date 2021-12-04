#pragma once

class Matrix;
class Vector;

#include <vector>
#include <stdexcept>

class Slice {
    public:
        Slice(const int& N, const int& M, const int& start, const int& step, const int& length) 
            : N(N), M(M), start(start), step(step), length(length) {}

        double dot(const Slice& s);
        double dot(const Vector& v);
        double norm();
        operator Vector() const;
        virtual const double& operator[](const int& j) const = 0;

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
        VectorSlice(const Vector* v);
        virtual ~VectorSlice() {}

        double& operator[](const int& i);
        const Vector* parent;
};

class MutableSlice : public Slice {
    public:
        MutableSlice(Matrix* const parent, const int& start, const int& step, const int& length);
        virtual ~MutableSlice() {}

        double& operator[](const int& i);
        const double& operator[](const int& i) const;
        MutableSlice& operator=(const Vector& v);
        MutableSlice& operator-=(const Slice& v);
        MutableSlice& operator-=(const Vector& v);
        MutableSlice& operator+=(const Slice& v);
        MutableSlice& operator+=(const Vector& v);

        Matrix* const parent;
};

class ConstSlice : public Slice {
    public:
        ConstSlice(const Matrix* parent, const int& start, const int& step, const int& length);
        virtual ~ConstSlice() {}

        const double& operator[](const int& i) const;

        const Matrix* const parent;
};

class ConstRow : public ConstSlice {
    public: 
        ConstRow(const Matrix* parent, const int& row);
        ConstRow(const ConstSlice& s) : ConstSlice(std::move(s)) {}
};

class ConstColumn : public ConstSlice {
    public: 
        ConstColumn(Matrix* const parent, const int& col);
        ConstColumn(const ConstSlice& s) : ConstSlice(std::move(s)) {}
};

class Row : public MutableSlice {
    public: 
        Row(Matrix* const parent, const int& row);
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
        Column(Matrix* const parent, const int& col);
        Column(MutableSlice& s) : MutableSlice(std::move(s)) {}

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