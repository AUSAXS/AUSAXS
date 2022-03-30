#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "Slice.h"

// A basic vector class. Sizes are checked before each operation, so an std::invalid_argument is thrown if they do not match.
template<typename T>
class Vector {
    public:
        /**
         * @brief Move constructor.
         */
        Vector(Vector<T>&& v) noexcept : _N(v._N), _data(std::move(v._data)) {}

        /**
         * @brief Copy constructor.
         */
        Vector(const Vector<T>& v) : _N(v.size()), _data(v._data) {}

        /**
         * @brief Construct a vector based on an initializer list.
         */
        Vector(const std::initializer_list<T> l) : _N(l.size()), _data(l) {}

        /**
         * @brief Construct a vector based on a std::vector. 
         */
        Vector(const std::vector<T>& v) : _N(v.size()), _data(v) {}

        /**
         * @brief Construct an empty vector of a given size. 
         */
        Vector(unsigned int n) : _N(n), _data(n) {}

        /**
         * @brief Default constructor.
         */
        Vector() : _N(0), _data(0) {}

        /**
         * @brief Destructor. 
         */
        virtual ~Vector() = default;

        // Assignment operator, w = v
        Vector<T>& operator=(const Vector<T>& v) {
            _N = v.N;
            _data = v.data;
            return *this;
        }

        // Slice assignment operator
        Vector<T>& operator=(const Slice<T>& s) {
			if (__builtin_expect(!(s.N == 1 || s.M == 1), false)) {throw std::invalid_argument("Only 1D slices can be assigned to vectors. Size: " + std::to_string(s.N) + ", " + std::to_string(s.M));}
			_N = std::max(s.N, s.M);
			_data = std::vector<double>(N);
			for (size_t i = 0; i < N; i++) {
				_data[i] = s[i];
			}
			return *this;
		}

        // Initializer list assignment operator
        Vector<T>& operator=(std::initializer_list<double> l) {
			_N = l.size();
			_data.assign(l);
			return *this;
		}

        // Plus operator, w + v
        template<typename Q>
        Vector<T> operator+(const Vector<Q>& v) const {
			compatibility_check(v);
			Vector w(N);
			std::transform(begin(), end(), v.begin(), w.begin(), std::plus<double>());
			return w;
		}

        // Minus operator, w - v
        template<typename Q>
        Vector<T> operator-(const Vector<Q>& v) const {
			compatibility_check(v);
			Vector w(N);
			std::transform(begin(), end(), v.begin(), w.begin(), std::minus<double>());
			return w;
		}

        // Minus operator, -w
        Vector<T> operator-() const {
			Vector w(N);
			std::transform(begin(), end(), w.begin(), std::negate<double>());
			return w;
		}

        // Vector multiplication, w*v
        template<typename Q>
        Vector<T> operator*(const Vector<Q>& v) const {
			compatibility_check(v);
			Vector w(N);
			std::transform(begin(), end(), v.begin(), w.begin(), std::multiplies<double>());
			return w;
		}

        // Scalar multiplication, w*a
        Vector<T> operator*(double a) const {
			Vector w(N);
			std::transform(begin(), end(), w.begin(), [&a](double x) {return x*a;});
			return w;
		}

        friend Vector<T> operator*(double a, const Vector<T>& v) {return v*a;}

        // Scalar division, w/a
        Vector<T> operator/(double a) const {
			Vector w(N);
			std::transform(begin(), end(), w.begin(), [&a](double x) {return x/a;});
			return w;
		}

        // Plus-assignment, w += v
        template<typename Q>
        Vector<T>& operator+=(const Vector<Q>& v) {
			compatibility_check(v);
			std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
			return *this;
		}

        // Minus-assignment, w -= v
        template<typename Q>
        Vector<T>& operator-=(const Vector<Q>& v) {
			compatibility_check(v);
			std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
			return *this;
		}

        // Scalar division-assignment, w /= a
        Vector<T>& operator/=(double a) {
			std::transform(begin(), end(), begin(), [&a](double x) {return x/a;});
			return *this;
		}

        // Scalar multiplication-assignment, w /= a
        Vector<T>& operator*=(double a) {
			std::transform(begin(), end(), begin(), [&a](double x) {return x*a;});
			return *this;
		}

        // Conversion to std::vector
        operator vector<T>() {
			return data;
		}

        // Read-only indexing, w[i]
        const T& operator[](unsigned int i) const {return data[i];}
        
        // Read/write indexing, w[i] = ...
        T& operator[](unsigned int i) {return _data[i];}

        // Approximate equality, w ~ v
        template<typename Q>
        bool operator==(const Vector<Q>& v) const {
			compatibility_check(v);
			Vector a = operator-(v); // difference vector
			return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
		}

        // Approximate inequality operator, w != v
        template<typename Q>
        bool operator!=(const Vector<Q>& v) const {return !operator==(v);}

        /**
         * @brief Get the dot product with another Vector.
         */
        template<typename Q>
        double dot(const Vector<Q>& v) const {
			compatibility_check(v);
			return std::inner_product(begin(), end(), v.begin(), 0.0);
		}

        /**
         * @brief Get the norm of this Vector.
         */
        double norm() const {return sqrt(dot(*this));}

        /**
         * @brief Get the magnitude of this Vector.
         */
        double magnitude() const {return norm();}

        /**
         * @brief Get the Euclidian distance to another Vector.
         */
        template<typename Q>
        double distance(const Vector<Q>& v) const {return sqrt(distance2(v));};

        /**
         * @brief Get the squared Euclidian distance to another Vector.
         */
        template<typename Q>
        double distance2(const Vector<Q>& v) const {
			compatibility_check(v);
			Vector w(N);
			std::transform(begin(), end(), v.begin(), w.begin(), [] (double x1, double x2) {return pow((x1-x2), 2);});
			return std::accumulate(w.begin(), w.end(), 0);
		}

        /**
         * @brief Get a copy of this Vector.
         */
        Vector copy() const {
			Vector w(N);
			std::copy(begin(), end(), w.begin());
			return w;
		}

        /**
         * @brief Get a string representation of this Vector.
         */
        std::string to_string(std::string message = "") const {
			if (message != "") {std::cout << message << std::endl;}
			std::stringstream s; s << "( ";
			for (const auto& e : data) {
				s << std::setprecision(8) << e << " ";
			}
			s << ")";
			return s.str();
		}

        /**
         * @brief Output the string representation of this Vector to a stream. 
         */
        friend std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {os << v.to_string(); return os;}

        // Read-only iterator
		const typename std::vector<T>::const_iterator begin() const {return data.begin();}

        // Read-only iterator
        const typename std::vector<T>::const_iterator end() const {return data.end();}

        // Read-write iterator
        typename std::vector<T>::iterator begin() {return _data.begin();}

        // Read-write iterator
        typename std::vector<T>::iterator end() {return _data.end();}

        /**
         * @brief Get the size of this Vector.
         */
        size_t size() const {return N;};

        /**
         * @brief Get the dimension of this Vector.
         */
        size_t dim() const {return size();}

        const size_t& N = _N; // read-only access to the dimension
        const std::vector<T>& data = _data; // read-only access to the data container

    protected:
        size_t _N;
        std::vector<T> _data;
        static constexpr double precision = 1e-9;

        // check if the vector is compatible with ours
        template<typename Q>
        void compatibility_check(const Vector<Q>& v) const {
			if (__builtin_expect(N != v.N, false)) {
				throw std::invalid_argument("Vector dimensions do not match (got: " + std::to_string(v.N) + ", expected: " + std::to_string(N) + ").");
			}
		}
};