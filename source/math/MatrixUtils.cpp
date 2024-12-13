/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <math/MatrixUtils.h>
#include <math/Vector3.h>
#include <math/Matrix.h>

#include <math.h>

using namespace ausaxs;

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