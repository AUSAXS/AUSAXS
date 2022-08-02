#pragma once

#include <initializer_list>

#include <math/Vector.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>

template<class T> 
class Vector3 : public Vector<T> {
	public:
		/**
		 * @brief Default constructor.
		 */
		Vector3() : Vector<T>(3) {}

		/**
		 * @brief Move constructor.
		 */
		Vector3(Vector3<T>&& v) noexcept : Vector<T>(std::move(v)) {} 

		/**
		 * @brief Copy constructor.
		 */
		Vector3(const Vector3<T>& v) : Vector<T>(v.data) {} 

		/**
		 * @brief List constructor.
		 */
		Vector3(const std::initializer_list<T> l) : Vector<T>(l) {}

		/**
		 * @brief Vector constructor.
		 */
		Vector3(const Vector<T>& v) : Vector<T>(v) {}

		/**
		 * @brief Constructor. 
		 * 
		 * @param x The x-coordinate.
		 * @param y The y-coordinate.
		 * @param z The z-coordinate. 
		 */
		Vector3(T x, T y, T z) : Vector<T>({x, y, z}) {}

		/**
		 * @brief Destructor.
		 */
		~Vector3() override = default;

		/**
		 * @brief Set this vector equal to another. 
		 */
		Vector3<T>& operator=(const Vector3<T>& v) {
			this->data = v.data;
			return *this;
		}

		/**
		 * @brief Set this vector equal to an initializer list. 
		 *        This allows the simple notation v = {1, 2, 3}.
		 */
		Vector3<T>& operator=(std::initializer_list<T> l) {
			this->N = l.size();
			this->data.assign(l);
			return *this;
		}

        // Plus operator, w + v
        template<typename Q>
        Vector3<T> operator+(const Vector3<Q>& v) const {return Vector3<T>(this->x() + v.x(), this->y() + v.y(), this->z() + v.z());}

        // Minus operator, w - v
        template<typename Q>
        Vector3<T> operator-(const Vector3<Q>& v) const {return Vector3<T>(this->x() - v.x(), this->y() - v.y(), this->z() - v.z());}

        // Minus operator, -w
        Vector3<T> operator-() const {return Vector3<T>(-this->x(), -this->y(), -this->z());}

        // Vector multiplication, w*v
        template<typename Q>
        Vector3<T> operator*(const Vector3<Q>& v) const {return Vector3<T>(this->x() * v.x(), this->y() * v.y(), this->z() * v.z());}

        // Scalar multiplication, w*a
        Vector3<T> operator*(double a) const {return Vector3<T>(this->x()*a, this->y()*a, this->z()*a);}

        friend Vector3<T> operator*(double a, const Vector3<T>& v) {return v*a;}

        // Scalar division, w/a
        Vector3<T> operator/(double a) const {return Vector3<T>(this->x()/a, this->y()/a, this->z()/a);}

        // Plus-assignment, w += v
        template<typename Q>
        Vector3<T>& operator+=(const Vector3<Q>& v) {
			this->x() += v.x();
			this->y() += v.y();
			this->z() += v.z();
			return *this;
		}

        // Minus-assignment, w -= v
        template<typename Q>
        Vector3<T>& operator-=(const Vector3<Q>& v) {
			this->x() -= v.x();
			this->y() -= v.y();
			this->z() -= v.z();
			return *this;
		}

        // Scalar division-assignment, w /= a
        Vector3<T>& operator/=(double a) {
			this->x() /= a;
			this->y() /= a;
			this->z() /= a;
			return *this;
		}

        // Scalar multiplication-assignment, w /= a
        Vector3<T>& operator*=(double a) {
			this->x() *= a;
			this->y() *= a;
			this->z() *= a;
			return *this;
		}

		/**
		 * @brief Get the Euclidian distance to another vector. 
		 */
		template<typename Q>
		double distance(const Vector3<Q>& v) const {return sqrt(pow(x()-v.x(), 2) + pow(y()-v.y(), 2) + pow(z()-v.z(), 2));}

		/**
		 * @brief Calculate the cross product of this vector with another. 
		 */
		template<typename Q>
		Vector3 cross(const Vector3<Q>& v) const {return {y()*v.z() - v.y()*z(), z()*v.x() - v.z()*x(), x()*v.y() - v.x()*y()};}

		/**
		 * @brief Rotate this vector around an axis. 
		 * 
		 * @param axis The rotation axis. 
		 * @param angle The angle to rotate. 
		 */
		template<typename Q>
		void rotate(const Vector3<Q>& axis, double angle) {
			Matrix R = matrix::rotation_matrix(axis, angle);
			rotate(R);
		}

		/**
		 * @brief Rotate this vector by a rotation matrix.
		 * 
		 * @param matrix The rotation matrix. 
		 */
		void rotate(const Matrix<double>& matrix) {
			*this = matrix*(*this);
		}

		/**
		 * @brief Output the string representation of this vector to a stream. 
		 */
		friend std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {os << v.to_string(); return os;}

		/**
		 * @brief Check if this vector is approximately equal to another. 
		 */
		template<typename Q>
		bool operator==(const Vector3<Q>& v) const {return abs(x()-v.x()) + abs(y()-v.y()) + abs(z()-v.z()) < this->precision;}

		/**
		 * @brief Check if this vector is compatible with another vector. Two 3D vectors are always compatible. 
		 */
		template<typename Q>
		void compatibility_check(const Vector3<Q>&) const {}

		/**
		 * @brief Normalize this vector to unit length.
		 */
		Vector3<double> normalize() {
			this->operator=(this->operator/(this->norm()));
			return *this;
		}

		/**
		 * @brief Get a copy of the normalized form of this vector.
		 */
		Vector3<double> normalize() const {
			return this->operator/(this->norm());
		}

		/**
		 * @brief Get a copy of the normalized form of this vector.
		 */
		Vector3<double> normalize_copy() const {
			return this->operator/(this->norm());
		}

		/**
		 * @brief Generate a complete 3D basis from this vector. 
		 *        Implementation based on the "frisvad" algorithm from https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf.
		 */
		std::tuple<Vector3<double>, Vector3<double>, Vector3<double>> generate_basis() {
			return vector3::generate_basis(*this);
		}

		inline double& x() {return this->data[0];}
		inline const double& x() const {return this->data[0];}

		inline double& y() {return this->data[1];}
		inline const double& y() const {return this->data[1];}

		inline double& z() {return this->data[2];}
		inline const double& z() const {return this->data[2];}
};