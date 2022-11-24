#pragma once

#include <utility/concepts.h>
#include <preprocessor.h>

#include <initializer_list>
#include <vector>
#include <string>
#include <iostream>

// A basic vector class. Sizes are checked before each operation, so an std::invalid_argument is thrown if they do not match.
template<numeric T>
class Vector {
    public:
        /**
         * @brief Move constructor.
         */
        Vector(Vector<T>&& v) noexcept : N(v.N), data(std::move(v.data)) {}

        /**
         * @brief Copy constructor.
         */
        Vector(const Vector<T>& v) : N(v.size()), data(v.data) {}

        /**
         * @brief Construct a vector based on an initializer list.
         */
        Vector(const std::initializer_list<T> l) : N(l.size()), data(l) {}

        /**
         * @brief Construct a vector based on a std::vector. 
         */
        Vector(const std::vector<T>& v) : N(v.size()), data(v) {}

        /**
         * @brief Construct an empty vector of a given size. 
         */
        Vector(unsigned int n) : N(n), data(n) {}

        /**
         * @brief Default constructor.
         */
        Vector() : N(0), data(0) {}

        /**
         * @brief Destructor. 
         */
        virtual ~Vector() = default;

        // Assignment operator, w = v
        Vector<T>& operator=(const Vector<T>& v);

        // Initializer list assignment operator
        Vector<T>& operator=(std::initializer_list<T> l);

        // Plus-assignment, w += v
        template<numeric Q>
        Vector<T>& operator+=(const Vector<Q>& v);

        // Minus-assignment, w -= v
        template<numeric Q>
        Vector<T>& operator-=(const Vector<Q>& v);

        // Scalar division-assignment, w /= a
        Vector<T>& operator/=(double a);

        // Scalar multiplication-assignment, w /= a
        Vector<T>& operator*=(double a);

        // Vector multiplication-assignment, w *= v
        template<numeric Q>
        Vector<T>& operator*=(const Vector<Q>& v);

        // Conversion to std::vector
        operator std::vector<T>();

        // Read-only indexing, w[i]
        const T& operator[](unsigned int i) const;
        
        // Read/write indexing, w[i] = ...
        T& operator[](unsigned int i);

        // Approximate equality, w ~ v
        template<numeric Q>
        bool operator==(const Vector<Q>& v) const;

        // Approximate inequality operator, w != v
        template<numeric Q>
        bool operator!=(const Vector<Q>& v) const;

        /**
         * @brief Get the dot product with another Vector.
         */
        template<numeric Q>
        double dot(const Vector<Q>& v) const;

        /**
         * @brief Get the norm of this Vector.
         */
        double norm() const;

        /**
         * @brief Get the magnitude of this Vector.
         */
        double magnitude() const;

        /**
         * @brief Get the Euclidian distance to another Vector.
         */
        template<numeric Q>
        double distance(const Vector<Q>& v) const;

        /**
         * @brief Get the squared Euclidian distance to another Vector.
         */
        template<numeric Q>
        double distance2(const Vector<Q>& v) const;

        /**
         * @brief Get a copy of this Vector.
         */
        Vector copy() const;

        /**
         * @brief Get a string representation of this Vector.
         */
        std::string to_string(std::string message = "") const;

        // Read-only iterator
		const typename std::vector<T>::const_iterator begin() const;

        // Read-only iterator
        const typename std::vector<T>::const_iterator end() const;

        // Read-write iterator
        typename std::vector<T>::iterator begin();

        // Read-write iterator
        typename std::vector<T>::iterator end();

        // Add data to the end of this Vector. 
        void push_back(T val);

        /**
         * @brief Get the size of this Vector.
         */
        size_t size() const;

        /**
         * @brief Get the dimension of this Vector.
         */
        size_t dim() const;

        void resize(unsigned int size);

        size_t N;
        std::vector<T> data;

    protected:
        static constexpr double precision = 1e-9;

        // check if the vector is compatible with ours
        template<numeric Q>
        void compatibility_check(const Vector<Q>& v) const;
};

template<numeric T, numeric Q>
Vector<T> operator+(Vector<T> left, const Vector<Q>& right) {return left += right;}

template<numeric T, numeric Q>
Vector<T> operator-(Vector<T> left, const Vector<Q>& right) {return left -= right;}

template<numeric T>
Vector<T> operator*(Vector<T> left, double right) {return left *= right;}

template<numeric T>
Vector<T> operator*(double left, Vector<T> right) {return right *= left;}

template<numeric T, numeric Q>
Vector<T> operator*(Vector<T> left, const Vector<Q>& right) {return left *= right;}

template<numeric T>
Vector<T> operator/(Vector<T> left, double right) {return left /= right;}

template<numeric T>
Vector<T> operator-(Vector<T> v) {return Vector<T>(v.size()) - v;}

template<numeric T> 
std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {os << v.to_string(); return os;}

#include <math/Vector.tpp>