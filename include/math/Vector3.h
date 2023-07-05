#pragma once

#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <math/Vector.h>
#include <utility/Concepts.h>
#include <utility/Exceptions.h>

#include <initializer_list>

template<numeric T> 
class Vector3 {
	public:
		/**
		 * @brief Default constructor.
		 */
		Vector3() : Vector3(0, 0, 0) {}

		/**
		 * @brief List constructor.
		 */
		Vector3(const std::initializer_list<T>& l) {
			*this = l;
		}

		Vector3(const Vector3<T>& v) : Vector3(v.copy()) {}
		Vector3(Vector3<T>&& v) : data(std::move(v.data)) {}

		Vector3(const Vector<T>& v) : Vector3(v.copy()) {}
		Vector3(Vector<T>&& v) {
			if (v.size() != 3) [[unlikely]] {
				throw except::invalid_argument("Vector3::Vector3: Vector must have size 3");
			}
			x() = v[0];
			y() = v[1];
			z() = v[2];
		}

		Vector3(const Matrix<T>& m) : Vector3(m.copy()) {}
		Vector3(Matrix<T>&& M) {
			if (M.N == 3 && M.M == 1) {
				x() = M.index(0, 0);
				y() = M.index(1, 0);
				z() = M.index(2, 0);
			} else if (M.N == 1 && M.M == 3) {
				x() = M.index(0, 0);
				y() = M.index(0, 1);
				z() = M.index(0, 2);
			} else {
				throw except::invalid_argument("Vector3::Vector3: Matrix must be 1x3 or 3x1");
			}
		}

		/**
		 * @brief Constructor. 
		 * 
		 * @param x The x-coordinate.
		 * @param y The y-coordinate.
		 * @param z The z-coordinate. 
		 */
		Vector3(T x, T y, T z) : data({x, y, z}) {}

		/**
		 * @brief Destructor.
		 */
		~Vector3() = default;

		T operator[](unsigned int i) const;

		T& operator[] (unsigned int i);

		/**
		 * @brief Set this vector equal to another. 
		 */
		Vector3<T>& operator=(const Vector3<T>& v);

		/**
		 * @brief Set this vector equal to an initializer list. 
		 *        This allows the simple notation v = {1, 2, 3}.
		 */
		Vector3<T>& operator=(std::initializer_list<T> l);

        // Plus-assignment, w += v
        template<numeric Q>
        Vector3<T>& operator+=(const Vector3<Q>& v);

        // Minus-assignment, w -= v
        template<numeric Q>
        Vector3<T>& operator-=(const Vector3<Q>& v);

        // Scalar division-assignment, w /= a
        Vector3<T>& operator/=(double a);

        // Scalar multiplication-assignment, w /= a
        Vector3<T>& operator*=(double a);

        // Vector multiplication-assignment, w /= v
		template<numeric Q>
        Vector3<T>& operator*=(const Vector3<Q>& v);

		template<numeric Q>
		bool operator==(const Vector3<Q>& v) const;

		template<numeric Q>
		bool operator!=(const Vector3<Q>& v) const;

		template<numeric Q>
		bool equals(const Vector3<Q>& v, double precision) const;

        /**
         * @brief Get the dot product with another Vector.
         */
        template<numeric Q>
        double dot(const Vector3<Q>& v) const;

        /**
         * @brief Get the norm of this Vector.
         */
        double norm() const;

        /**
         * @brief Get the magnitude of this Vector.
         */
        double magnitude() const;

		/**
		 * @brief Get the Euclidian distance to another vector. 
		 */
		template<numeric Q>
		double distance(const Vector3<Q>& v) const;

        /**
         * @brief Get the squared Euclidian distance to another Vector.
         */
        template<numeric Q>
        double distance2(const Vector3<Q>& v) const;

		/**
		 * @brief Calculate the cross product of this vector with another. 
		 */
		template<numeric Q>
		Vector3 cross(const Vector3<Q>& v) const;

		/**
		 * @brief Rotate this vector around an axis. 
		 * 
		 * @param axis The rotation axis. 
		 * @param angle The angle to rotate. 
		 */
		template<numeric Q>
		void rotate(const Vector3<Q>& axis, double angle);

		/**
		 * @brief Rotate this vector by a rotation matrix.
		 * 
		 * @param matrix The rotation matrix. 
		 */
		void rotate(const Matrix<double>& matrix);

		/**
		 * @brief Normalize this vector to unit length.
		 */
		Vector3<T>& normalize();

		/**
		 * @brief Generate a complete 3D basis from this vector. 
		 *        Implementation based on the "frisvad" algorithm from https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf.
		 */
		std::tuple<Vector3<double>, Vector3<double>, Vector3<double>> generate_basis();

        /**
         * @brief Get a string representation of this Vector.
         */
        std::string to_string(std::string message = "") const;

		// Conversion to std::vector
        operator std::vector<T>();

        // Conversion to Vector
        operator Vector<T>();

        // Conversion to Matrix
        operator Matrix<T>();

		T& x();
		const T& x() const;

		T& y();
		const T& y() const;

		T& z();
		const T& z() const;

		size_t size() const;

		Vector3<T> copy() const;

		typename std::array<T, 3>::iterator begin();
		typename std::array<T, 3>::iterator end();
		const typename std::array<T, 3>::const_iterator begin() const;
		const typename std::array<T, 3>::const_iterator end() const;

        static constexpr double precision = 1e-9;

	private:
		std::array<T, 3> data;
};

template<numeric T, numeric Q>
Vector3<Q> operator*(const Matrix<T>& M, const Vector3<Q>& v) {
    #if (SAFE_MATH)
        if (M.M != v.size()) [[unlikely]] {
            throw except::invalid_argument("Vector3::operator*: Invalid matrix dimensions (got: " + std::to_string(M.M) + ", expected: " + std::to_string(v.size()) + "]).");
        }
    #endif

    Vector3<Q> w;
    for (size_t row = 0; row < v.size(); ++row) {
        for (size_t col = 0; col < M.M; ++col) {
			w[row] += M.index(row, col) * v[col];
        }
    }
    return w;
}

template<numeric T, numeric Q>
bool operator==(const Vector3<T>& v, const Vector<Q>& w) {
	#if (SAFE_MATH)
		if (v.size() != w.size()) [[unlikely]] {
			throw except::invalid_argument("Vector3::operator*: Invalid vector dimensions (got: " + std::to_string(v.size()) + ", expected: " + std::to_string(w.size()) + "]).");
		}
	#endif

	return abs(v.x() - w[0]) + abs(v.y() - w[1]) + abs(v.z() - w[2]) < Vector3<T>::precision;
}

template<numeric T, numeric Q>
bool operator==(const Vector<T>& v, const Vector3<Q>& w) {return w == v;}

template<numeric T, numeric Q>
bool operator!=(const Vector3<T>& v, const Vector<Q>& w) {return !(v == w);}

template<numeric T, numeric Q>
bool operator!=(const Vector<T>& v, const Vector3<Q>& w) {return !(w == v);}

template<numeric T, numeric Q>
Vector3<T> operator+(Vector3<T> left, const Vector3<Q>& right) {return left += right;}

template<numeric T, numeric Q>
Vector3<T> operator-(Vector3<T> left, const Vector3<Q>& right) {return left -= right;}

template<numeric T>
Vector3<T> operator-(Vector3<T> v) {return Vector3<T>() - v;}

template<numeric T>
Vector3<T> operator*(Vector3<T> left, double right) {return left *= right;}

template<numeric T>
Vector3<T> operator*(double left, Vector3<T> right) {return right *= left;}

template<numeric T, numeric Q>
Vector3<T> operator*(Vector3<T> left, const Vector3<Q>& right) {return left *= right;}

template<numeric T>
Vector3<T> operator/(Vector3<T> left, double right) {return left /= right;}

template<numeric T>
Vector3<T> operator/(double left, Vector3<T> right) {return Vector3<T>(left/right.x(), left/right.y(), left/right.z());}

template<numeric T>
std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {os << v.to_string(); return os;}

#include <math/Vector3.tpp>