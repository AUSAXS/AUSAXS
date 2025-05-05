#pragma once

#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <math/Vector.h>
#include <math/MathConcepts.h>
#include <math/MathTypeTraits.h>

#include <initializer_list>
#include <stdexcept>
#include <array>

namespace ausaxs {
	template<numeric T> 
	class Vector3 final {
		public:
			Vector3() = default;
			Vector3(const Vector3<T>& v) = default;
			Vector3(Vector3<T>&& v) = default;
			Vector3& operator=(const Vector3<T>& v) = default;
			Vector3& operator=(Vector3<T>&& v) = default;

			Vector3(const std::initializer_list<T>& l) {
				*this = l;
			}

			template<typename Q> requires std::is_floating_point_v<Q>
			Vector3(const Vector3<Q>& v) {
				data[0] = static_cast<T>(v[0]);
				data[1] = static_cast<T>(v[1]);
				data[2] = static_cast<T>(v[2]);
			}

			Vector3(const Vector<T>& v) {
				assert(v.size() == 3 && "Vector3: Vector must have size 3");
				std::copy_n(v.begin(), 3, data.begin());
			}

			Vector3(Vector<T>&& v) {
				assert(v.size() == 3 && "Vector3: Vector must have size 3");
				std::copy_n(v.begin(), 3, data.begin());
			}

			Vector3(const Matrix<T>& M) {
				assert(M.data.size() == 3 && "Vector3: Matrix must have size 3");
				std::copy_n(M.begin(), 3, data.begin());
			}

			Vector3(Matrix<T>&& M) {
				assert(M.data.size() == 3 && "Vector3: Matrix must have size 3");
				std::copy_n(M.begin(), 3, data.begin());
			}

			/**
			 * @brief Constructor. 
			 * 
			 * @param x The x-coordinate.
			 * @param y The y-coordinate.
			 * @param z The z-coordinate. 
			 */
			Vector3(T x, T y, T z) : data({x, y, z}) {}

			T operator[](unsigned int i) const;

			T& operator[] (unsigned int i);

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

			operator std::vector<T>();
			operator Vector<T>();
			operator Matrix<T>();

			T& x();
			T& y();
			T& z();
			const T& x() const;
			const T& y() const;
			const T& z() const;

			size_t size() const;

			Vector3<T> copy() const;

			template<size_t i> T& get();
			template<size_t i> const T& get() const;

			typename std::array<T, 3>::iterator begin();
			typename std::array<T, 3>::iterator end();
			const typename std::array<T, 3>::const_iterator begin() const;
			const typename std::array<T, 3>::const_iterator end() const;

			static constexpr double precision = 1e-6;

		private:
			std::array<T, 3> data;
	};

	template<numeric T, numeric Q>
	Vector3<Q> operator*(const Matrix<T>& M, const Vector3<Q>& v) {
		#if (SAFE_MATH)
			if (M.M != v.size()) [[unlikely]] {
				throw std::invalid_argument("Vector3::operator*: Invalid matrix dimensions (got: " + std::to_string(M.M) + ", expected: " + std::to_string(v.size()) + "]).");
			}
		#endif

		return {
			M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
			M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
			M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]
		};
	}

	template<numeric T, numeric Q>
	bool operator==(const Vector3<T>& v, const Vector<Q>& w) {
		#if (SAFE_MATH)
			if (v.size() != w.size()) [[unlikely]] {
				throw std::invalid_argument("Vector3::operator*: Invalid vector dimensions (got: " + std::to_string(v.size()) + ", expected: " + std::to_string(w.size()) + "]).");
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

	template<numeric T>
	Vector3<T> operator+(Vector3<T> left, double right) {
		Vector3<T> w = std::move(left);
		w.x() += right;
		w.y() += right;
		w.z() += right;
		return w;
	}

	template<numeric T, numeric Q>
	Vector3<T> operator-(Vector3<T> left, const Vector3<Q>& right) {return left -= right;}

	template<numeric T>
	Vector3<T> operator-(Vector3<T> left, double right) {
		Vector3<T> w = std::move(left);
		w.x() -= right;
		w.y() -= right;
		w.z() -= right;
		return w;
	}

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
}

namespace std {
	template<> struct tuple_size<ausaxs::Vector3<double>> : std::integral_constant<size_t, 3> {};
	template<> struct tuple_size<ausaxs::Vector3<float>> : std::integral_constant<size_t, 3> {};
	template<> struct tuple_size<ausaxs::Vector3<int>> : std::integral_constant<size_t, 3> {};
	template<> struct tuple_element<0, ausaxs::Vector3<double>> {using type = double;};	
	template<> struct tuple_element<1, ausaxs::Vector3<double>> {using type = double;};
	template<> struct tuple_element<2, ausaxs::Vector3<double>> {using type = double;};
	template<> struct tuple_element<0, ausaxs::Vector3<float>> {using type = float;};
	template<> struct tuple_element<1, ausaxs::Vector3<float>> {using type = float;};
	template<> struct tuple_element<2, ausaxs::Vector3<float>> {using type = float;};
	template<> struct tuple_element<0, ausaxs::Vector3<int>> {using type = int;};
	template<> struct tuple_element<1, ausaxs::Vector3<int>> {using type = int;};
	template<> struct tuple_element<2, ausaxs::Vector3<int>> {using type = int;};
}

#include <math/Vector3.tpp>
static_assert(sizeof(ausaxs::Vector3<double>) == 24, "Vector3 size is not 24 bytes");
static_assert(std::is_trivial_v<ausaxs::Vector3<double>>, "Vector3 is not trivial");
static_assert(std::is_standard_layout_v<ausaxs::Vector3<double>>, "Vector3 is not standard layout");
static_assert(supports_nothrow_move_v<ausaxs::Vector3<double>>, "Vector3 should support nothrow move semantics.");