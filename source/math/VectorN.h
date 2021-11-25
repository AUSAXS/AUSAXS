#include <initializer_list>
#include <vector>

template <typename T>
class VectorN {
public:
    VectorN(std::initializer_list<T> l) : data(l), N(l.size()) {}
    VectorN(vector<T> v) : data(v), N(v.size()) {}

    VectorN operator+(const VectorN& w) const;
    VectorN& operator+=(const VectorN& w) const;
    VectorN operator-(const VectorN& w) const;
    VectorN& operator-=(const VectorN& w) const;
    VectorN operator*(const double& a) const;
    VectorN operator/(const double& a) const;
    double dot(const VectorN& w) const;
    double norm() const;
    void distance(const VectorN& w) const;

private:
    const int N;
    std::vector<T> data;
};

template <typename T>
VectorN<T> VectorN<T>::operator+(const VectorN<T>& w) const {
    std::vector<T> v(N);
    for (int i = 0; i < N, i++) {v[i] = data[i] + w[i];}
    return v;
}

template <typename T>
VectorN<T>& VectorN<T>::operator+=(const VectorN<T>& w) const {
    for (int i = 0; i < N; i++) {data[i] += w[i];}
    return *this;
}

template <typename T>
VectorN<T> VectorN<T>::operator-(const VectorN<T>& w) const {
    std::vector<T> v(N);
    for (int i = 0; i < N, i++) {v[i] = data[i] - w[i];}
    return v;
}

template <typename T>
VectorN<T>& VectorN<T>::operator-=(const VectorN<T>& w) const {
    for (int i = 0; i < N; i++) {data[i] -= w[i];}
    return *this;
}

template <typename T>
VectorN<T> VectorN<T>::operator*(const double& a) const {
    std::vector<T> v(N);
    for (int i = 0; i < N; i++) {v[i] = data[i]*a;}
    return v;
}

template <typename T>
VectorN<T> VectorN<T>::operator/(const double& a) const {
    std::vector<T> v(N);
    for (int i = 0; i < N; i++) {v[i] = data[i]/a;}
    return v;
}

template <typename T>
double VectorN<T>::dot(const VectorN<T>& w) const {
    double a = 0;
    for (int i = 0; i < N; i++) {a += data[i]*w[i];}
    return a;
}

template <typename T>
double VectorN<T>::norm() const {return dot(*this);}