#pragma once

#include <math/MathConcepts.h>
#include <math/MathTypeTraits.h>

#include <initializer_list>
#include <vector>
#include <string>
#include <iostream>

namespace ausaxs {
    // A basic vector class. Sizes are checked before each operation, so an std::invalid_argument is thrown if they do not match.
    template<numeric T>
    class Vector final {
        public:
            Vector() = default;
            Vector(const Vector<T>& v) = default;
            Vector(Vector<T>&& v) = default;
            Vector& operator=(const Vector<T>& v) = default;
            Vector& operator=(Vector<T>&& v) = default;

            /**
             * @brief Construct a vector based on an initializer list.
             */
            Vector(const std::initializer_list<T> l) : data(l) {}

            /**
             * @brief Construct a vector based on a std::vector. 
             */
            Vector(const std::vector<T>& v) : data(v) {}

            /**
             * @brief Construct a vector based on a std::vector. 
             */
            Vector(std::vector<T>&& v) : data(std::move(v)) {}

            /**
             * @brief Construct an empty vector of a given size. 
             */
            Vector(unsigned int n) : data(n) {}

            Vector<T>& operator=(std::initializer_list<T> l);

            template<numeric Q>
            Vector<T>& operator+=(const Vector<Q>& v);

            template<numeric Q>
            Vector<T>& operator-=(const Vector<Q>& v);

            Vector<T>& operator/=(double a);
            Vector<T>& operator*=(double a);

            template<numeric Q>
            Vector<T>& operator*=(const Vector<Q>& v);

            // Conversion to std::vector. This is a O(N) operation.
            operator std::vector<T>();

            const T& operator[](unsigned int i) const;
            T& operator[](unsigned int i);

            // Approximate equality, w ~ v
            template<numeric Q>
            bool operator==(const Vector<Q>& v) const;

            // Approximate inequality operator, w !~ v
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
             * @brief Get a string representation of this Vector.
             */
            std::string to_string() const;

            const typename std::vector<T>::const_iterator begin() const;
            const typename std::vector<T>::const_iterator end() const;
            typename std::vector<T>::iterator begin();
            typename std::vector<T>::iterator end();

            // Add data to the end of this Vector. 
            void push_back(T val);

            /**
             * @brief Get the size of this Vector.
             */
            unsigned int size() const;

            /**
             * @brief Get the dimension of this Vector.
             */
            unsigned int dim() const;

            /**
             * @brief Resize this Vector to the given size.
             */
            void resize(unsigned int size);

            /**
             * @brief Check if this Vector is empty.
             */
            bool empty() const;

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
}

#include <math/Vector.tpp>
static_assert(std::is_standard_layout_v<ausaxs::Vector<double>>, "Vector is not standard layout");
static_assert(supports_nothrow_move_v<ausaxs::Vector<double>>, "Vector should support nothrow move semantics.");