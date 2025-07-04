// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/slices/Slice.h>
#include <math/Vector.h>
#include <math/MathConcepts.h>
#include <math/MathTypeTraits.h>

#include <initializer_list>

namespace ausaxs {
    /**
     * @brief A simple matrix class. 
     *        This uses a single 1D array to efficiently store the data, and offers a variety of vector operations.
     *        Note that this class is _not_ optimized for speed! If efficiency is a concern, consider using e.g. Eigen. 
     */
    template<numeric Q>
    class Matrix final {
        public: 
            Matrix() : N(0), M(0) {}
            Matrix(const Matrix<Q>& A) = default;
            Matrix(Matrix<Q>&& A) = default;
            Matrix& operator=(const Matrix<Q>& A) = default;
            Matrix& operator=(Matrix<Q>&& A) = default;

            template<typename T> requires std::is_floating_point_v<T>
            Matrix(const Matrix<T>& A) : N(A.N), M(A.M), data(N*M) {
                std::copy(A.begin(), A.end(), data.begin());
            }

            /**
             * @brief Construct a Matrix based on a nested initializer list. The lists must be of the same size. 
             */
            Matrix(std::initializer_list<std::initializer_list<Q>> l);

            /**
             * @brief Construct a Matrix based on a list of column vectors. The vectors must be of the same size.
             */
            Matrix(const std::vector<std::vector<Q>>& cols);

            /**
             * @brief Construct a Matrix based on a vector.
             */
            Matrix(const Vector<Q>& v);

            /**
             * @brief Construct an empty Matrix of a given size. 
             */
            Matrix(unsigned int n, unsigned int m);

            /**
             * @brief Add a new row at the end of the matrix.
             */
            void push_back(const std::vector<double>& r);

            /**
             * @brief Get the identity matrix of a given dimension. 
             */
            static Matrix<Q> identity(unsigned int dim);

            Matrix<Q> operator-() const;
            Matrix<Q>& operator*=(double a);
            Matrix<Q>& operator/=(double a);
            template<numeric R> Matrix<Q>& operator+=(const Matrix<R>& A);
            template<numeric R> Matrix<Q>& operator-=(const Matrix<R>& A);

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

            const ConstRow<Q> operator[](unsigned int i) const;
            MutableRow<Q> operator[](unsigned int i);
            const ConstColumn<Q> col(unsigned int j) const;
            MutableColumn<Q> col(unsigned int j);
            const ConstRow<Q> row(unsigned int i) const;
            MutableRow<Q> row(unsigned int i);

            // Approximate equality operator
            template<numeric R>
            bool operator==(const Matrix<R>& A) const;

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

            const Q& operator()(unsigned int i, unsigned int j) const;
            Q& operator()(unsigned int i, unsigned int j);

            const Q& index(unsigned int i, unsigned int j) const;
            Q& index(unsigned int i, unsigned int j);

            const typename std::vector<Q>::const_iterator begin() const;
            const typename std::vector<Q>::const_iterator end() const;

            typename std::vector<Q>::iterator begin();
            typename std::vector<Q>::iterator end();

            std::string to_string() const;

            unsigned int N, M;
            std::vector<Q> data;
            static constexpr double precision = 1e-9;

        private: 
            /**
             * @brief Check if the matrix is compatible with ours.
             *        This check can be disabled by setting the macro SAFE_MATH to 0.
             */
            template<numeric R>
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

    template<numeric Q, numeric R>
    Matrix<Q> operator+(Matrix<Q> left, const Matrix<R>& right) {return left += right;}

    template<numeric Q, numeric R>
    Matrix<Q> operator-(Matrix<Q> left, const Matrix<R>& right) {return left -= right;} 

    template<numeric Q>
    Matrix<Q> operator*(Matrix<Q> left, double right) {return left *= right;}

    template<numeric Q>
    Matrix<Q> operator*(double left, Matrix<Q> right) {return right *= left;}

    template<numeric Q>
    Matrix<Q> operator/(Matrix<Q> left, double right) {return left /= right;}

    template<numeric Q, numeric R>
    Vector<Q> operator*(const Matrix<Q>& A, const Vector<R>& v) {
        #if (SAFE_MATH)
            if (A.M != v.size()) [[unlikely]] {
                throw std::invalid_argument("Matrix::operator*: Invalid matrix dimensions (got: " + std::to_string(v.size()) + ", expected: " + std::to_string(A.M) + "]).");
            }
        #endif

        Vector<Q> w(A.N);
        for (unsigned int row = 0; row < A.N; ++row) {
            Q sum = Q();
            for (unsigned int col = 0; col < A.M; ++col) {
                sum += v[col]*A[row][col];
            }
            w[row] = sum;
        }
        return w;
    }

    template<numeric Q, numeric R>
    Matrix<Q> operator*(const Matrix<Q>& A, const Matrix<R>& B) {
        #if (SAFE_MATH)
            if (A.M != B.N) [[unlikely]] {
                throw std::invalid_argument("Matrix::operator*: Invalid matrix dimensions (got: " + std::to_string(A.M) + ", " + std::to_string(A.N) + ", expected: " + std::to_string(B.N) + ", " + std::to_string(B.M) + "]).");
            }
        #endif

        Matrix<Q> C(A.N, B.M);
        for (unsigned int row = 0; row < C.N; row++) {
            for (unsigned int col = 0; col < C.M; col++) {
                Q sum = Q();
                for (unsigned int inner = 0; inner < A.M; inner++) {
                    sum += A[row][inner]*B[inner][col];
                }
                C[row][col] = sum;
            }
        }
        return C;
    }
}

#include <math/Matrix.tpp>
static_assert(supports_nothrow_move_v<ausaxs::Matrix<double>>, "Matrix should support nothrow move semantics.");