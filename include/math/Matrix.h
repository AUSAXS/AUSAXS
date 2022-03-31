#pragma once

#include <iostream>
#include <iterator>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iomanip>

#include <math/Slice.h>
#include <math/Vector.h>
#include <math/Vector3.h>
#include <math/LUPDecomposition.h>

template<typename Q>
class Matrix {
    public: 
        /**
         * @brief Move constructor. 
         */
        Matrix(Matrix<Q>&& A) noexcept : N(A.N), M(A.M), data(std::move(A.data)) {}

        /**
         * @brief Copy constructor.
         */
        Matrix(const Matrix<Q>& A) : N(A.N), M(A.M), data(A.data) {}

        /**
         * @brief Construct a Matrix based on a nested initializer list. The lists must be of the same size. 
         */
        Matrix(std::initializer_list<std::initializer_list<Q>> l) : N(l.size()), M(l.begin()->size()) {
            for (const auto& row : l) {
                if (__builtin_expect(row.size() != M, false)) {throw std::invalid_argument("Malformed matrix: columns must be of equal size!");}
                for (const auto& e : row) {
                    data.push_back(e);
                }
            }
        }

        /**
         * @brief Construct a Matrix based on a vector.
         */
        Matrix(const Vector<Q>& v) : N(v.N), M(1), data(v.data) {}

        /**
         * @brief Construct an empty Matrix of a given size. 
         */
        Matrix(int n, int m) : N(n), M(m), data(N*M) {} 

        /**
         * @brief Default constructor.
         */
        Matrix() : N(0), M(0), data(0) {}

        /**
         * @brief Destructor.
         */
        ~Matrix() = default;

        /**
         * @brief Get the identity matrix of a given dimension. 
         */
        static Matrix<Q> identity(unsigned int dim) {
            Matrix<Q> A(dim, dim);
            for (unsigned int i = 0; i < dim; i++) {
                A[i][i] = 1;
            }
            return A;
        } 

        // Assignment operator, B = A
        Matrix<Q>& operator=(const Matrix<Q>& A) {
            N = A.N; M = A.M;
            data = A.data;
            return *this;
        }

        // Plus operator, B + A
        template<typename R>
        Matrix<Q> operator+(const Matrix<R>& A) const {
            compatibility_check(A);
            Matrix<Q> B(N, M);
            std::transform(begin(), end(), A.begin(), B.begin(), std::plus<double>());
            return B;
        }

        // Minus operator, B - A
        template<typename R>
        Matrix<Q> operator-(const Matrix<R>& A) const {
            compatibility_check(A);
            Matrix<Q> B(N, M);
            std::transform(begin(), end(), A.begin(), B.begin(), std::minus<double>());
            return B;
        }

        // Negation operator, -A
        Matrix<Q> operator-() const {
            Matrix<Q> A(N, M);
            std::transform(begin(), end(), A.begin(), std::negate<double>());
            return A;
        }

        // Scalar multiplication, B*a
        Matrix<Q> operator*(double a) const {
            Matrix<Q> A(N, M);
            std::transform(begin(), end(), A.begin(), [&a] (const double& e) {return e*a;});
            return A;
        }

        // Scalar division, B/a
        Matrix<Q> operator/(double a) const {
            Matrix<Q> A(N, M);
            std::transform(begin(), end(), A.begin(), [&a] (const double& e) {return e/a;});
            return A;
        }

        // Plus-assignment, B += A
        template<typename R>
        Matrix<Q>& operator+=(const Matrix<R>& A) {
            compatibility_check(A);
            std::transform(begin(), end(), A.begin(), begin(), std::plus<double>());
            return *this;
        }

        // Minus-assignment, B -= A
        template<typename R>
        Matrix<Q>& operator-=(const Matrix<R>& A) {
            compatibility_check(A);
            std::transform(begin(), end(), A.begin(), begin(), std::minus<double>());
            return *this;
        }

        // Vector multiplication, A*v
        template<typename R>
        friend Vector<Q> operator*(const Matrix<Q>& A, const Vector<R>& v) {
            if (__builtin_expect(A.M != v.N, false)) {
                throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(v.N) + ", expected: " + std::to_string(A.M) + "]).");
            }
            Vector<Q> w(A.N);
            for (size_t row = 0; row < A.N; ++row) {
                for (size_t col = 0; col < A.M; ++col) {
                    w[row] += v[col]*A[row][col];
                }
            }
            return w;
        }

        // Matrix multiplication, A*B
        template<typename R>
        friend Matrix<Q> operator*(const Matrix<Q>& A, const Matrix<R>& B) {
            if (__builtin_expect(A.M != B.N, false)) {
                throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(A.M) + ", " + std::to_string(A.N) + 
                    ", expected: " + std::to_string(B.N) + ", " + std::to_string(B.M) + "]).");
            }
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

        // Read-only indexer
        const ConstRow<Q> operator[](unsigned int i) const {return row(i);}
        // Read-write indexer
        Row<Q> operator[](unsigned int i) {return row(i);}

        // Read-only column indexer
        const ConstColumn<Q> col(unsigned int j) const {return ConstColumn<Q>(data, N, M, j);}
        // Read-write column indexer
        Column<Q> col(unsigned int j) {return Column<Q>(data, N, M, j);}

        // Read-only row indexer
        const ConstRow<Q> row(unsigned int i) const {return ConstRow<Q>(data, N, M, i);}
        // Read-write row indexer
        Row<Q> row(unsigned int i) {return Row<Q>(data, N, M, i);}

        // Approximate equality, B ~ A
        template<typename R>
        bool operator==(const Matrix<R>& A) const {
            compatibility_check(A);
            Matrix diff = operator-(A); // difference matrix
            return std::accumulate(diff.begin(), diff.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
        }

        // Approximate inequality operator, w != v
        template<typename R>
        bool operator!=(const Matrix<R>& A) const {return !operator==(A);}

        /**
         * @brief Get the determinant of this Matrix.
         */
        double det() const {
            if (__builtin_expect(N != M, false)) {throw std::invalid_argument("Error in matrix determinant: Matrix is not square.");}
            LUPDecomposition decomp(*this);
            return decomp.determinant();
        }

        /**
         * @brief Copy this Matrix. 
         */
        Matrix<Q> copy() const {
            Matrix A(N, M);
            A.data.assign(data.begin(), data.end());
            return A;
        }

        /**
         * @brief Get the transpose of this Matrix.
         */
        Matrix<Q> T() const {
            Matrix A(M, N);
            for (size_t row = 0; row < A.N; ++row) {
                for (size_t col = 0; col < A.M; ++col) {
                    A[row][col] = index(col, row);
                }
            }
            return A;
        }

        /**
         * @brief Get the transpose of this Matrix.
         */
        Matrix<Q> transpose() const {
            return T();
        }

        // Read-only indexer
        const Q& index(unsigned int i, unsigned int j) const {return data[M*i + j];}

        // Read-write indexer
        Q& index(unsigned int i, unsigned int j) {return data[M*i + j];}

        // Read-only iterator
        const typename std::vector<Q>::const_iterator begin() const {return data.begin();}

        // Read-only iterator
        const typename std::vector<Q>::const_iterator end() const {return data.end();}

        // Read-write iterator
        typename std::vector<Q>::iterator begin() {return data.begin();}

        // Read-write iterator
        typename std::vector<Q>::iterator end() {return data.end();}

        /**
         * @brief Format and print this Matrix to the terminal.
         */
        void print(std::string message = "") const {
            if (!message.empty()) {std::cout << message << std::endl;}
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
        static Matrix<double> rotation_matrix(double alpha, double beta, double gamma) {
            double cosa = cos(alpha), cosb = cos(beta), cosg = cos(gamma);
            double sina = sin(alpha), sinb = sin(beta), sing = sin(gamma);
            double sinasinb = sina*sinb, cosasinb = cosa*sinb;

            return Matrix{{cosb*cosg, sinasinb*cosg - cosa*sing, cosasinb*cosg + sina*sing}, 
                            {cosb*sing, sinasinb*sing + cosa*cosg, cosasinb*sing - sina*cosg},
                            {-sinb,     sina*cosb,                 cosa*cosb}};
        }

        /**
         * @brief Generate a 3x3 rotation matrix from a rotation axis and an angle around this axis. 
         *        This uses the Euler-Rodrigues formulation.
         * @param axis The rotation axis.
         * @param angle The rotation angle.
         */
        static Matrix<double> rotation_matrix(const Vector3& axis, double angle) {
            Vector3 ax = axis.normalize_copy();
            double a = cos(angle/2);
            double b = sin(angle/2);
            double c = b;
            double d = b;
            b *= ax.x;
            c *= ax.y;
            d *= ax.z;

            double aa = a*a, bb = b*b, cc = c*c, dd = d*d;
            double bc = b*c, ad = a*d, ac = a*c, ab = a*b, bd = b*d, cd = c*d;

            Matrix R{{aa+bb-cc-dd, 2*(bc-ad),   2*(bd+ac)}, 
                    {2*(bc+ad),   aa+cc-bb-dd, 2*(cd-ab)},
                    {2*(bd-ac),   2*(cd+ab),   aa+dd-bb-cc}};
            return R;
        }

        size_t N, M;
        std::vector<Q> data;
        static constexpr double precision = 1e-9;

    private: 
        // check if the matrix is compatible with ours
        template<typename R>
        void compatibility_check(const Matrix<R>& A) const {
            if (__builtin_expect(N != A.N || M != A.M, false)) {
                throw std::invalid_argument("Matrix dimensions do not match (got: [" + std::to_string(N) + ", " + std::to_string(M) + "] and [" + 
                    std::to_string(A.N) + ", " + std::to_string(A.M) + "]).");
            }
        }
};