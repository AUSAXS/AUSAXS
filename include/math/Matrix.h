#pragma once

#include <iostream>
#include <iterator>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iomanip>

#include <math/slices/ConstSlice.h>
#include <math/slices/MutableSlice.h>
#include <math/Vector.h>
#include <math/LUPDecomposition.h>

#define SAFE_MATH true

#include <signal.h>
template<typename Q> // using Q to avoid conflict with transpose method T()
class Matrix {
    public: 
        /**
         * @brief Move constructor. 
         */
        Matrix(Matrix<Q>&& A) noexcept;

        /**
         * @brief Copy constructor.
         */
        Matrix(const Matrix<Q>& A);

        /**
         * @brief Construct a Matrix based on a nested initializer list. The lists must be of the same size. 
         */
        Matrix(std::initializer_list<std::initializer_list<Q>> l);

        /**
         * @brief Construct a Matrix based on a series of vectors. The vectors must be of the same size. 
         */
        Matrix(std::vector<std::vector<Q>> v);

        /**
         * @brief Construct a Matrix based on a vector.
         */
        Matrix(const Vector<Q>& v);

        /**
         * @brief Construct an empty Matrix of a given size. 
         */
        Matrix(int n, int m);

        /**
         * @brief Default constructor.
         */
        Matrix();

        /**
         * @brief Destructor.
         */
        virtual ~Matrix();

        /**
         * @brief Add a new row at the end of the matrix.
         */
        void push_back(std::vector<double> r);

        /**
         * @brief Get the identity matrix of a given dimension. 
         */
        static Matrix<Q> identity(unsigned int dim);

        Matrix<Q>& operator=(const Matrix<Q>& A);
        Matrix<Q> operator-() const;
        Matrix<Q>& operator*=(double a);
        Matrix<Q>& operator/=(double a);
        template<typename R> Matrix<Q>& operator+=(const Matrix<R>& A);
        template<typename R> Matrix<Q>& operator-=(const Matrix<R>& A);

        /**
         * @brief Extend the number of rows by the specified amount. Data in the old rows is preserved. 
         *        Complexity: O(N*M)
         */
        void extend(int n);

        /**
         * @brief Resize the matrix to a new shape. 
         *        If the number of columns is changed, the content of the new matrix is undefined. 
         *        Otherwise the new rows are filled with zeros, leaving the old rows intact. 
         *        Complexity: O(N*M)
         */
        void resize(int n, int m);

        // Read-only indexer
        const ConstRow<Q> operator[](unsigned int i) const;

        // Read-write indexer
        Row<Q> operator[](unsigned int i);

        // Read-only column indexer
        const ConstColumn<Q> col(unsigned int j) const;

        // Read-write column indexer
        Column<Q> col(unsigned int j);

        // Read-only row indexer
        const ConstRow<Q> row(unsigned int i) const;

        // Read-write row indexer
        Row<Q> row(unsigned int i);

        // Approximate equality operator
        template<typename R>
        bool operator==(const Matrix<R>& A) const;

        // Approximate inequality operator
        template<typename R>
        bool operator!=(const Matrix<R>& A) const;

        /**
         * @brief Get the determinant of this Matrix.
         */
        double det() const;

        /**
         * @brief Copy this Matrix. 
         */
        Matrix<Q> copy() const;

        /**
         * @brief Get the transpose of this Matrix.
         */
        Matrix<Q> T() const;

        /**
         * @brief Get the transpose of this Matrix.
         */
        Matrix<Q> transpose() const;

        // Read-only indexer
        const Q& index(unsigned int i, unsigned int j) const;

        // Read-write indexer
        Q& index(unsigned int i, unsigned int j);

        // Read-only iterator
        const typename std::vector<Q>::const_iterator begin() const;

        // Read-only iterator
        const typename std::vector<Q>::const_iterator end() const;

        // Read-write iterator
        typename std::vector<Q>::iterator begin();

        // Read-write iterator
        typename std::vector<Q>::iterator end();

        /**
         * @brief Format and print this Matrix to the terminal.
         */
        void print(std::string message = "") const;

        size_t N, M;
        std::vector<Q> data;
        static constexpr double precision = 1e-9;

    private: 
        /**
         * @brief Check if the matrix is compatible with ours.
         *        This check can be disabled by setting the macro SAFE_MATH to 0.
         */
        template<typename R>
        void compatibility_check(const Matrix<R>& A) const;

        /**
         * @brief Check if the number of columns is compatible with ours. 
         *        This check can be disabled by setting the macro SAFE_MATH to 0.
         */
        void compatibility_check_N(unsigned int N) const;

        /**
         * @brief Check if the number of rows is compatible with ours. 
         *        This check can be disabled by setting the macro SAFE_MATH to 0.
         */
        void compatibility_check_M(unsigned int M) const;
};

template<typename Q, typename R>
Matrix<Q> operator+(Matrix<Q> left, const Matrix<R>& right) {return left += right;}

template<typename Q, typename R>
Matrix<Q> operator-(Matrix<Q> left, const Matrix<R>& right) {return left -= right;} 

template<typename Q>
Matrix<Q> operator*(Matrix<Q> left, double right) {return left *= right;}

template<typename Q>
Matrix<Q> operator*(double left, Matrix<Q> right) {return right *= left;}

template<typename Q>
Matrix<Q> operator/(Matrix<Q> left, double right) {return left /= right;}

template<typename Q, typename R>
Vector<Q> operator*(const Matrix<Q>& A, const Vector<R>& v) {
    #if (SAFE_MATH)
        if (__builtin_expect(A.M != v.N, false)) {
            throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(v.N) + ", expected: " + std::to_string(A.M) + "]).");
        }
    #endif

    Vector<Q> w(A.N);
    for (size_t row = 0; row < A.N; ++row) {
        for (size_t col = 0; col < A.M; ++col) {
            w[row] += v[col]*A[row][col];
        }
    }
    return w;
}

template<typename Q, typename R>
Matrix<Q> operator*(const Matrix<Q>& A, const Matrix<R>& B) {
    #if (SAFE_MATH)
        if (__builtin_expect(A.M != B.N, false)) {
            throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(A.M) + ", " + std::to_string(A.N) + ", expected: " + std::to_string(B.N) + ", " + std::to_string(B.M) + "]).");
        }
    #endif

    Matrix<Q> C(A.N, B.M);
    for (size_t row = 0; row < C.N; row++) {
        for (size_t col = 0; col < C.M; col++) {
            for (size_t inner = 0; inner < A.M; inner++) {
                C[row][col] += A[row][inner]*B[inner][col];
            }
        }
    }
    return C;
}

#include <math/Matrix.tpp>