#pragma once

#include <math/Vector.h>
#include <utility/Exceptions.h>

#include <functional>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <math.h>

template<numeric T>
Vector<T>& Vector<T>::operator=(const Vector<T>& v) {
    N = v.N;
    data = v.data;
    return *this;
}

template<numeric T>
Vector<T>& Vector<T>::operator=(std::initializer_list<T> l) {
    N = l.size();
    data.assign(l);
    return *this;
}

template<numeric T> template<numeric Q>
Vector<T>& Vector<T>::operator+=(const Vector<Q>& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
    return *this;
}

template<numeric T> template<numeric Q>
Vector<T>& Vector<T>::operator-=(const Vector<Q>& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
    return *this;
}

template<numeric T>
Vector<T>& Vector<T>::operator/=(double a) {
    std::transform(begin(), end(), begin(), [&a] (T x) {return x/a;});
    return *this;
}

template<numeric T>
Vector<T>& Vector<T>::operator*=(double a) {
    std::transform(begin(), end(), begin(), [&a] (T x) {return x*a;});
    return *this;
}

template<numeric T> template<numeric Q>
Vector<T>& Vector<T>::operator*=(const Vector<Q>& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::multiplies<T>());
    return *this;
}

template<numeric T>
Vector<T>::operator std::vector<T>() {
    return std::move(data);
}

template<numeric T>
const T& Vector<T>::operator[](unsigned int i) const {return data[i];}

template<numeric T>
T& Vector<T>::operator[](unsigned int i) {return data[i];}

template<numeric T> template<numeric Q>
bool Vector<T>::operator==(const Vector<Q>& v) const {
    compatibility_check(v);
    Vector<T> a = *this - v; // difference vector
    return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, T x) {return sum + abs(x);}) < precision;
}

template<numeric T> template<numeric Q>
bool Vector<T>::operator!=(const Vector<Q>& v) const {return !operator==(v);}

template<numeric T> template<numeric Q>
double Vector<T>::dot(const Vector<Q>& v) const {
    compatibility_check(v);
    return std::inner_product(begin(), end(), v.begin(), 0.0);
}

template<numeric T>
double Vector<T>::norm() const {return sqrt(dot(*this));}

template<numeric T>
double Vector<T>::magnitude() const {return norm();}

template<numeric T> template<numeric Q>
double Vector<T>::distance(const Vector<Q>& v) const {return sqrt(distance2(v));}

template<numeric T> template<numeric Q>
double Vector<T>::distance2(const Vector<Q>& v) const {
    compatibility_check(v);
    Vector<T> w(N);
    std::transform(begin(), end(), v.begin(), w.begin(), [] (T x1, Q x2) {return pow((x1-x2), 2);});
    return std::accumulate(w.begin(), w.end(), 0);
}

template<numeric T>
Vector<T> Vector<T>::copy() const {
    Vector w(N);
    std::copy(begin(), end(), w.begin());
    return w;
}

template<numeric T>
std::string Vector<T>::to_string() const {
    std::stringstream s; s << "( ";
    for (const auto& e : data) {
        s << std::setprecision(8) << e << " ";
    }
    s << ")";
    return s.str();
}

template<numeric T>
const typename std::vector<T>::const_iterator Vector<T>::begin() const {return data.cbegin();}

template<numeric T>
const typename std::vector<T>::const_iterator Vector<T>::end() const {return data.cend();}

template<numeric T>
typename std::vector<T>::iterator Vector<T>::begin() {return data.begin();}

template<numeric T>
typename std::vector<T>::iterator Vector<T>::end() {return data.end();}

template<numeric T>
void Vector<T>::push_back(T val) {data.push_back(val); N++;}

template<numeric T>
unsigned Vector<T>::size() const {return N;}

template<numeric T>
unsigned Vector<T>::dim() const {return size();}

template<numeric T>
void Vector<T>::resize(unsigned int size) {
    N = size;
    data.resize(size);
}

template<numeric T> template<numeric Q>
void Vector<T>::compatibility_check(const Vector<Q>& v) const {
    #if (SAFE_MATH)
        if (N != v.N) [[unlikely]] {
            throw except::invalid_argument("Vector::compatibility_check: Vector dimensions do not match (got: " + std::to_string(v.N) + ", expected: " + std::to_string(N) + ").");
        }
    #endif
}