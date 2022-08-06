#include <math/MatrixUtils.h>
#include <math/Vector3.h>
#include <math/Matrix.h>

Matrix<double> matrix::rotation_matrix(double alpha, double beta, double gamma) {
    double cosa = cos(alpha), cosb = cos(beta), cosg = cos(gamma);
    double sina = sin(alpha), sinb = sin(beta), sing = sin(gamma);
    double sinasinb = sina*sinb, cosasinb = cosa*sinb;

    return Matrix{{cosb*cosg, sinasinb*cosg - cosa*sing, cosasinb*cosg + sina*sing}, 
                    {cosb*sing, sinasinb*sing + cosa*cosg, cosasinb*sing - sina*cosg},
                    {-sinb,     sina*cosb,                 cosa*cosb}};
}

Matrix<double> matrix::rotation_matrix(const Vector3<double>& axis, double angle) {
    auto ax = axis; 
    ax.normalize();
    double a = cos(angle/2);
    double b = sin(angle/2);
    double c = b;
    double d = b;
    b *= ax.x();
    c *= ax.y();
    d *= ax.z();

    double aa = a*a, bb = b*b, cc = c*c, dd = d*d;
    double bc = b*c, ad = a*d, ac = a*c, ab = a*b, bd = b*d, cd = c*d;

    Matrix R{{aa+bb-cc-dd, 2*(bc-ad),   2*(bd+ac)}, 
            {2*(bc+ad),   aa+cc-bb-dd, 2*(cd-ab)},
            {2*(bd-ac),   2*(cd+ab),   aa+dd-bb-cc}};
    return R;
}

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