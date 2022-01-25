#pragma once

class Vector;

#include <vector>
#include <stdexcept>

using std::vector;

/**
 * @brief \class Slice
 * 
 * A Slice is a type of cut of a Matrix, typically either in the form of a single row or column. 
 * The idea with this class is to give mutable access to such a structure. 
 */
class Slice {
  public:

    /**
     * @brief Constructor. 
     * 
     * @param N The number of rows of this Slice. 
     * @param M The number of columns of this Slice. 
     * @param start The start location in the raw data array. 
     * @param step The number of elements to skip between each index in the raw data array. 
     * @param length The total number of elements this Slice can access. 
     */
    Slice(const int N, const int M, const int start, const int step, const int length) 
        : N(N), M(M), start(start), step(step), length(length) {}

    /**
     * @brief Get the dot product with another Slice. 
     */
    double dot(const Slice& s);

    /**
     * @brief Get the dot product with a Vector.
     */
    double dot(const Vector& v);

    /**
     * @brief Get the norm of this Slice. 
     */
    double norm();

    /**
     * @brief Cast this Slice into a Vector. 
     */
    operator Vector() const;

    /**
     * @brief Mutable indexer in this Slice.
     * 
     * @param j The index to access.
     */
    virtual const double& operator[](const int j) const = 0;

    /**
     * @brief Format and print the contents of this Slice to the terminal. 
     */
    void print() const;

    const size_t N, M;
    const int start, step, length;
};

// Add a Vector to this Slice. 
Vector operator+(const Slice& s, const Vector& v);

// Subtract a Vector from this Slice.
Vector operator-(const Slice& s, const Vector& v);

// Multiply each element of this Slice with a constant. 
Vector operator*(const Slice& s, const double a);

// Divide each element of this Slice with a constant.
Vector operator/(const Slice& s, const double a);

// Check if this Slice is equal to a given Vector.
bool operator==(const Slice& s, const Vector& v);

// Check if this Slice is not equal to a given Vector.
bool operator!=(const Slice& s, const Vector& v);

// VectorSlice class - unused. 
class VectorSlice : public Slice {
  public: 
    VectorSlice(const vector<double>& data);
    virtual ~VectorSlice() {}

    double& operator[](const int i);
    const vector<double>& data;
};

/**
 * @brief \class MutableSlice. 
 * 
 * This class represents a mutable slice. 
 */
class MutableSlice : public Slice {
  public:
    MutableSlice(vector<double>& data, const int N, const int M, const int start, const int step, const int length);
    virtual ~MutableSlice() {}

    double& operator[](const int i);
    const double& operator[](const int i) const;
    MutableSlice& operator=(const Vector& v);
    MutableSlice& operator-=(const Slice& v);
    MutableSlice& operator-=(const Vector& v);
    MutableSlice& operator+=(const Slice& v);
    MutableSlice& operator+=(const Vector& v);

    vector<double>& data;
};

/**
 * @brief \class ConstSlice.
 * 
 * This class represents read-only slices.
 */
class ConstSlice : public Slice {
  public:
    ConstSlice(const vector<double>& data, const int N, const int M, const int start, const int step, const int length);
    virtual ~ConstSlice() {}

    const double& operator[](const int i) const;

    const vector<double>& data;
};

class ConstRow : public ConstSlice {
  public: 
    ConstRow(const vector<double>& data, const int N, const int M, const int row);
    ConstRow(const ConstSlice& s) : ConstSlice(std::move(s)) {}
};

class ConstColumn : public ConstSlice {
  public: 
    ConstColumn(const vector<double>& data, const int N, const int M, const int col);
    ConstColumn(const ConstSlice& s) : ConstSlice(std::move(s)) {}
};

class Row : public MutableSlice {
  public: 
    Row(vector<double>& data, const int N, const int M, const int row);
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
    Column(vector<double>& data, const int N, const int M, const int col);
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