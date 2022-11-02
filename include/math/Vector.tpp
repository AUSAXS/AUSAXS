#pragma once

#include <math/Vector.h>

template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& v) {
    N = v.N;
    data = v.data;
    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator=(std::initializer_list<T> l) {
    N = l.size();
    data.assign(l);
    return *this;
}

template<typename T> template<typename Q>
Vector<T>& Vector<T>::operator+=(const Vector<Q>& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
    return *this;
}

template<typename T> template<typename Q>
Vector<T>& Vector<T>::operator-=(const Vector<Q>& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator/=(double a) {
    std::transform(begin(), end(), begin(), [&a] (T x) {return x/a;});
    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator*=(double a) {
    std::transform(begin(), end(), begin(), [&a] (T x) {return x*a;});
    return *this;
}

template<typename T> template<typename Q>
Vector<T>& Vector<T>::operator*=(const Vector<Q>& v) {
    compatibility_check(v);
    std::transform(begin(), end(), v.begin(), begin(), std::multiplies<T>());
    return *this;
}

template<typename T>
Vector<T>::operator std::vector<T>() {
    return data;
}

template<typename T>
const T& Vector<T>::operator[](unsigned int i) const {return data[i];}

template<typename T>
T& Vector<T>::operator[](unsigned int i) {return data[i];}

template<typename T> template<typename Q>
bool Vector<T>::operator==(const Vector<Q>& v) const {
    compatibility_check(v);
    Vector<T> a = *this - v; // difference vector
    return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, T x) {return sum + abs(x);}) < precision;
}

template<typename T> template<typename Q>
bool Vector<T>::operator!=(const Vector<Q>& v) const {return !operator==(v);}

template<typename T> template<typename Q>
double Vector<T>::dot(const Vector<Q>& v) const {
    compatibility_check(v);
    return std::inner_product(begin(), end(), v.begin(), 0.0);
}

template<typename T>
double Vector<T>::norm() const {return sqrt(dot(*this));}

template<typename T>
double Vector<T>::magnitude() const {return norm();}

template<typename T> template<typename Q>
double Vector<T>::distance(const Vector<Q>& v) const {return sqrt(distance2(v));}

template<typename T> template<typename Q>
double Vector<T>::distance2(const Vector<Q>& v) const {
    compatibility_check(v);
    Vector<T> w(N);
    std::transform(begin(), end(), v.begin(), w.begin(), [] (T x1, Q x2) {return pow((x1-x2), 2);});
    return std::accumulate(w.begin(), w.end(), 0);
}

template<typename T>
Vector<T> Vector<T>::copy() const {
    Vector w(N);
    std::copy(begin(), end(), w.begin());
    return w;
}

template<typename T>
std::string Vector<T>::to_string(std::string message) const {
    if (message != "") {std::cout << message << std::endl;}
    std::stringstream s; s << "( ";
    for (const auto& e : data) {
        s << std::setprecision(8) << e << " ";
    }
    s << ")";
    return s.str();
}

template<typename T>
const typename std::vector<T>::const_iterator Vector<T>::begin() const {return data.begin();}

template<typename T>
const typename std::vector<T>::const_iterator Vector<T>::end() const {return data.end();}

template<typename T>
typename std::vector<T>::iterator Vector<T>::begin() {return data.begin();}

template<typename T>
typename std::vector<T>::iterator Vector<T>::end() {return data.end();}

template<typename T>
void Vector<T>::push_back(T val) {data.push_back(val); N++;}

template<typename T>
size_t Vector<T>::size() const {return N;}

template<typename T>
size_t Vector<T>::dim() const {return size();}

template<typename T>
void Vector<T>::resize(unsigned int size) {
    N = size;
    data.resize(size);
}

template<typename T> template<typename Q>
void Vector<T>::compatibility_check(const Vector<Q>& v) const {
    #if (SAFE_MATH)
        if (__builtin_expect(N != v.N, false)) {
            throw std::invalid_argument("Vector dimensions do not match (got: " + std::to_string(v.N) + ", expected: " + std::to_string(N) + ").");
        }
    #endif
}