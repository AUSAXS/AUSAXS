#pragma once

#include <vector>
#include "Vector.h"

class Matrix;

template <typename T>
class Slice {
    public:
        Slice(T parent) : parent(parent) {}

        virtual std::vector<double> copy() const;
        explicit virtual operator Vector() const;
        double dot(const Vector& v) {return Vector().dot(v);}
        double norm() {return Vector().norm();}
        Vector operator+(const Vector& v) const {return Vector().operator+(v);}
        Vector operator-(const Vector& v) const {return Vector().operator-(v);}
        Vector operator*(const double& a) const {return Vector().operator*(a);}
        Vector operator/(const double& a) const {return Vector().operator/(a);}
        bool operator==(const Vector& v) const {return Vector().operator==(v);}
        bool operator!=(const Vector& v) const {return !operator==(v);}

        int index;
        const int N, M;
        T parent;
};

class RowSlice : public Slice<Matrix* const> {
    public:
        RowSlice(Matrix* const parent, const int& row);

        RowSlice& operator=(const RowSlice& s);
        RowSlice& operator=(const Vector& v);
        std::vector<double> copy() const;
        const double& operator[](const int& j) const;
        double& operator[](const int& j);
        operator Vector() const;

    private: 
        const int row;
};

class ColumnSlice : public Slice<Matrix* const> {
    public:
        ColumnSlice(Matrix* const parent, const int& col);

        ColumnSlice& operator=(const ColumnSlice& s);
        ColumnSlice& operator=(const Vector& v);
        std::vector<double> copy() const;
        const double& operator[](const int& j) const;
        double& operator[](const int& j);
        operator Vector() const;

    private: 
        const int col;
};

class ConstRowSlice : public Slice<const Matrix*> {
    public:
        ConstRowSlice(const Matrix* parent, const int& row);

        std::vector<double> copy() const;
        const double& operator[](const int& j) const;
        operator Vector() const;

    private: 
        const int row;
};

class ConstColumnSlice : public Slice<const Matrix*> {
    public:
        ConstColumnSlice(Matrix* const parent, const int& col);

        std::vector<double> copy() const;
        const double& operator[](const int& j) const;
        operator Vector() const;

    private: 
        const int col;
};