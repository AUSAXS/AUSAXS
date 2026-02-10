// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <math/MatrixUtils.h>
#include <math/Vector3.h>
#include <math/Matrix.h>

#include <math.h>
#include <cassert>

using namespace ausaxs;

template<numeric T>
bool is_rotation_matrix(const Matrix<T>& R) {
    Matrix<T> should_be_identity = R.transpose() * R;
    return std::accumulate(should_be_identity.begin(), should_be_identity.end(), 0.0, [] (double acc, double val) {return acc + std::abs(val);}) - 3 < 1e-6;
}

template<numeric T>
Vector3<T> matrix::euler_angles(const Matrix<T>& R) {
    assert(is_rotation_matrix(R) && "Input matrix is not a valid rotation matrix.");
    T c = std::hypot(R.index(0, 0), R.index(1, 0));
    T beta = std::atan2(-R.index(2, 0), c);
    T alpha, gamma;
    if (1e-9 < c) {
        alpha = std::atan2(R(2, 1), R(2, 2));
        gamma = std::atan2(R(1, 0), R(0, 0));
    } else { // gimbal lock
        alpha = 0.0;
        gamma = std::atan2(-R(0, 1), R(1, 1));
    }
    return {alpha, beta, gamma};
}

template<numeric T>
Matrix<T> matrix::rotation_matrix(T alpha, T beta, T gamma) {
    double cosa = std::cos(alpha), cosb = std::cos(beta), cosg = std::cos(gamma);
    double sina = std::sin(alpha), sinb = std::sin(beta), sing = std::sin(gamma);
    double sinasinb = sina*sinb, cosasinb = cosa*sinb;

    return Matrix{
        {static_cast<T>(cosb*cosg), static_cast<T>(sinasinb*cosg - cosa*sing), static_cast<T>(cosasinb*cosg + sina*sing)}, 
        {static_cast<T>(cosb*sing), static_cast<T>(sinasinb*sing + cosa*cosg), static_cast<T>(cosasinb*sing - sina*cosg)},
        {static_cast<T>(-sinb),     static_cast<T>(sina*cosb),                 static_cast<T>(cosa*cosb)}
    };
}

template<numeric T>
Matrix<T> matrix::rotation_matrix(const Vector3<T>& angles) {return rotation_matrix(angles.x(), angles.y(), angles.z());}

template<numeric T>
Matrix<T> matrix::rotation_matrix(const Vector3<T>& axis, double angle) {
    auto ax = axis; 
    ax.normalize();
    double a = std::cos(angle/2);
    double b = std::sin(angle/2);
    double c = b;
    double d = b;
    b *= ax.x();
    c *= ax.y();
    d *= ax.z();

    double aa = a*a, bb = b*b, cc = c*c, dd = d*d;
    double bc = b*c, ad = a*d, ac = a*c, ab = a*b, bd = b*d, cd = c*d;

    Matrix R{
        {static_cast<T>(aa+bb-cc-dd), static_cast<T>(2*(bc-ad)),   static_cast<T>(2*(bd+ac))}, 
        {static_cast<T>(2*(bc+ad)),   static_cast<T>(aa+cc-bb-dd), static_cast<T>(2*(cd-ab))},
        {static_cast<T>(2*(bd-ac)),   static_cast<T>(2*(cd+ab)),   static_cast<T>(aa+dd-bb-cc)}
    };
    return R;
}

template Matrix<double> matrix::rotation_matrix(double alpha, double beta, double gamma);
template Matrix<double> matrix::rotation_matrix(const Vector3<double>& angles);
template Matrix<double> matrix::rotation_matrix(const Vector3<double>& axis, double angle);
template Matrix<float> matrix::rotation_matrix(float alpha, float beta, float gamma);
template Matrix<float> matrix::rotation_matrix(const Vector3<float>& angles);
template Matrix<float> matrix::rotation_matrix(const Vector3<float>& axis, double angle);
template Vector3<double> matrix::euler_angles(const Matrix<double>& R);
template Vector3<float> matrix::euler_angles(const Matrix<float>& R);

Matrix<double> matrix::identity(unsigned int dim) {
    Matrix<double> A(dim, dim);
    for (unsigned int i = 0; i < dim; i++) {
        A.index(i, i) = 1;
    }
    return A;
} 

std::tuple<Vector3<double>, Vector3<double>, Vector3<double>> vector3::generate_basis(const Vector3<double>& v) {
    Vector3 n = v;
    n.normalize();

    // Handle the singularity
    if (n.z() < -0.9999999) { 
        Vector3<double> b1(0, -1, 0);
        Vector3<double> b2(-1, 0, 0);
        return std::make_tuple(n, b1, b2);
    }
    const float a = 1/(1 + n.z());
    const float b = -n.x()*n.y()*a;
    Vector3<double> b1(1-n.x()*n.x()*a, b, -n.x());
    Vector3<double> b2(b, 1-n.y() * n.y()*a, -n.y());
    return std::make_tuple(n, b1, b2);
}