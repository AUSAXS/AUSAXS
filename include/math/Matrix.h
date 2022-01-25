#pragma once

#include <iostream>
#include <iterator>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iomanip>

#include "Slice.h"
#include "Vector.h"

class Matrix {
    friend class Row;
    friend class Column;
    friend class MutableSlice;

  public: 
    /**
     * @brief Move constructor. 
     */
    Matrix(Matrix&& A) noexcept : _N(A.N), _M(A.M), _data(std::move(A._data)) {}

    /**
     * @brief Copy constructor.
     */
    Matrix(const Matrix& A) : _N(A.N), _M(A.M), _data(A.data) {}

    /**
     * @brief Construct a Matrix based on a nested initializer list. The lists must be of the same size. 
     */
    Matrix(std::initializer_list<std::initializer_list<double>> l) : _N(l.size()), _M(l.begin()->size()) {
        for (const auto& row : l) {
            if (__builtin_expect(row.size() != M, false)) {throw std::invalid_argument("Malformed matrix: columns must be of equal size!");}
            for (const auto& e : row) {
                _data.push_back(e);
            }
        }
    }

    /**
     * @brief Construct a Matrix based on a vector.
     */
    Matrix(const Vector& v);

    /**
     * @brief Construct an empty Matrix of a given size. 
     */
    Matrix(const int n, const int m) : _N(n), _M(m), _data(N*M) {} 

    /**
     * @brief Default constructor.
     */
    Matrix() : _N(0), _M(0), _data(0) {}

    /**
     * @brief Destructor.
     */
    ~Matrix() = default;

    /**
     * @brief Get the identity matrix of a given dimension. 
     */
    static Matrix identity(const int dim) {
        Matrix A(dim, dim);
        for (int i = 0; i < dim; i++) {
            A[i][i] = 1;
        }
        return A;
    } 

    // Assignment operator, B = A
    Matrix& operator=(const Matrix& A) {
        _N = A.N; _M = A.M;
        _data = A.data;
        return *this;
    }

    // Plus operator, B + A
    Matrix operator+(const Matrix& A) const {
        compatibility_check(A);
        Matrix B(N, M);
        std::transform(begin(), end(), A.begin(), B.begin(), std::plus<double>());
        return B;
    }

    // Minus operator, B - A
    Matrix operator-(const Matrix& A) const {
        compatibility_check(A);
        Matrix B(N, M);
        std::transform(begin(), end(), A.begin(), B.begin(), std::minus<double>());
        return B;
    }

    // Negation operator, -A
    Matrix operator-() const {
        Matrix A(N, M);
        std::transform(begin(), end(), A.begin(), std::negate<double>());
        return A;
    }

    // Scalar multiplication, B*a
    Matrix operator*(const double a) const {
        Matrix A(N, M);
        std::transform(begin(), end(), A.begin(), [&a] (const double& e) {return e*a;});
        return A;
    }

    // Scalar division, B/a
    Matrix operator/(const double a) const {
        Matrix A(N, M);
        std::transform(begin(), end(), A.begin(), [&a] (const double& e) {return e/a;});
        return A;
    }

    // Plus-assignment, B += A
    Matrix& operator+=(const Matrix& A) {
        compatibility_check(A);
        std::transform(begin(), end(), A.begin(), begin(), std::plus<double>());
        return *this;
    }

    // Minus-assignment, B -= A
    Matrix& operator-=(const Matrix& A) {
        compatibility_check(A);
        std::transform(begin(), end(), A.begin(), begin(), std::minus<double>());
        return *this;
    }

    // Vector multiplication, A*v
    friend Vector operator*(const Matrix& A, const Vector& v) {
        if (__builtin_expect(A.M != v.N, false)) {
            throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(v.N) + ", expected: " + std::to_string(A.M) + "]).");
        }
        Vector w(A.N);
        for (size_t row = 0; row < A.N; ++row) {
            for (size_t col = 0; col < A.M; ++col) {
                w[row] += v[col]*A[row][col];
            }
        }
        return w;
    }

    // Matrix multiplication, A*B
    friend Matrix operator*(const Matrix& A, const Matrix& B) {
        if (__builtin_expect(A.M != B.N, false)) {
            throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(A.M) + ", " + std::to_string(A.N) + 
                ", expected: " + std::to_string(B.N) + ", " + std::to_string(B.M) + "]).");
        }
        Matrix C(A.N, B.M);
        for (size_t row = 0; row < C.N; row++) {
            for (size_t col = 0; col < C.M; col++) {
                for (size_t inner = 0; inner < A.M; inner++) {
                    C[row][col] += A[row][inner]*B[inner][col];
                }
            }
        }
        return C;
    }

    // Read-only indexer
    const ConstRow operator[](const int i) const;

    // Read-only column indexer
    const ConstColumn col(const int j) const;

    // Read-only row indexer
    const ConstRow row(const int i) const;
    
    // Read-write indexer
    Row operator[](const int i);

    // Read-write row indexer
    Row row(const int i);

    // Read-write column indexer
    Column col(const int j);

    // Approximate equality, B ~ A
    bool operator==(const Matrix& A) const {
        compatibility_check(A);
        Matrix diff = operator-(A); // difference matrix
        return std::accumulate(diff.begin(), diff.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
    }

    // Approximate inequality operator, w != v
    bool operator!=(const Matrix& A) const {return !operator==(A);}

    /**
     * @brief Get the determinant of this Matrix.
     */
    double det() const;

    /**
     * @brief Copy this Matrix. 
     */
    Matrix copy() const {
        Matrix A(N, M);
        A._data.assign(data.begin(), data.end());
        return A;
    }

    /**
     * @brief Get the transpose of this Matrix.
     */
    Matrix T() const {
        Matrix A(M, N);
        for (size_t row = 0; row < A.N; ++row) {
            for (size_t col = 0; col < A.M; ++col) {
                A[row][col] = index(col, row);
            }
        }
        return A;
    }

    // Read-only indexer
    const double& index(const int i, const int j) const {return data[M*i + j];}

    // Read-write indexer
    double& index(const int i, const int j) {return _data[M*i + j];}

    // Read-only iterator
    const std::vector<double>::const_iterator begin() const {return data.begin();}

    // Read-only iterator
    const std::vector<double>::const_iterator end() const {return data.end();}

    // Read-write iterator
    std::vector<double>::iterator begin() {return _data.begin();}

    // Read-write iterator
    std::vector<double>::iterator end() {return _data.end();}

    /**
     * @brief Format and print this Matrix to the terminal.
     */
    void print(const std::string& message = "") const {
        if (message != "") {std::cout << message << std::endl;}
        for (size_t i = 0; i < N; i++) {
            std::cout << "\t" << std::setprecision(3);
            for (size_t j = 0; j < M; j++) {
                std::cout << std::setw(8) << index(i, j);
            }
            std::cout << std::endl;
        }
    }

    /**
     * @brief Generate a 3x3 extrinsic rotation matrix.
     */
    static Matrix rotation_matrix(const double alpha, const double beta, const double gamma) {
        double cosa = cos(alpha), cosb = cos(beta), cosg = cos(gamma);
        double sina = sin(alpha), sinb = sin(beta), sing = sin(gamma);
        double sinasinb = sina*sinb, cosasinb = cosa*sinb;

        return Matrix{{cosb*cosg, sinasinb*cosg - cosa*sing, cosasinb*cosg + sina*sing}, 
                      {cosb*sing, sinasinb*sing + cosa*cosg, cosasinb*sing - sina*cosg},
                      {-sinb,     sina*cosb,                 cosa*cosb}};
    }

    const size_t &N = _N, &M = _M; // read-only access to the dimensions
    const std::vector<double>& data = _data; // read-only access to the data container

  private: 
    size_t _N, _M;
    std::vector<double> _data;
    static constexpr double precision = 1e-9;

    // check if the matrix is compatible with ours
    virtual void compatibility_check(const Matrix& A) const {
        if (__builtin_expect(N != A.N || M != A.M, false)) {
            throw std::invalid_argument("Matrix dimensions do not match (got: [" + std::to_string(N) + ", " + std::to_string(M) + "] and [" + 
                std::to_string(A.N) + ", " + std::to_string(A.M) + "]).");
        }
    }
};