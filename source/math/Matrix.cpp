#include "math/Matrix.h"
#include "math/Vector.h"
#include "math/Vector3.h"
#include "math/LUPDecomposition.h"

Matrix::Matrix(const Vector& v) : _N(v.N), _M(1), _data(v.data) {} // vector --> matrix constructor

double Matrix::det() const {
    if (__builtin_expect(N != M, false)) {throw std::invalid_argument("Error in matrix determinant: Matrix is not square.");}
    LUPDecomposition decomp(*this);
    return decomp.determinant();
}

Matrix Matrix::rotation_matrix(const Vector3& axis, double angle) {
    Vector3 ax = axis.normalize_copy();
    double a = cos(angle/2);
    double b = sin(angle/2);
    double c = b;
    double d = b;
    b *= ax.x;
    c *= ax.y;
    d *= ax.z;

    double aa = a*a, bb = b*b, cc = c*c, dd = d*d;
    double bc = b*c, ad = a*d, ac = a*c, ab = a*b, bd = b*d, cd = c*d;

    Matrix R{{aa+bb-cc-dd, 2*(bc-ad),   2*(bd+ac)}, 
            {2*(bc+ad),   aa+cc-bb-dd, 2*(cd-ab)},
            {2*(bd-ac),   2*(cd+ab),   aa+dd-bb-cc}};
    return R;
}

Matrix Matrix::rotation_matrix(double alpha, double beta, double gamma) {
    double cosa = cos(alpha), cosb = cos(beta), cosg = cos(gamma);
    double sina = sin(alpha), sinb = sin(beta), sing = sin(gamma);
    double sinasinb = sina*sinb, cosasinb = cosa*sinb;

    return Matrix{{cosb*cosg, sinasinb*cosg - cosa*sing, cosasinb*cosg + sina*sing}, 
                    {cosb*sing, sinasinb*sing + cosa*cosg, cosasinb*sing - sina*cosg},
                    {-sinb,     sina*cosb,                 cosa*cosb}};
}

void Matrix::print(std::string message) const {
    if (message != "") {std::cout << message << std::endl;}
    for (size_t i = 0; i < N; i++) {
        std::cout << "\t" << std::setprecision(3);
        for (size_t j = 0; j < M; j++) {
            std::cout << std::setw(8) << index(i, j);
        }
        std::cout << std::endl;
    }
}

Matrix Matrix::T() const {
    Matrix A(M, N);
    for (size_t row = 0; row < A.N; ++row) {
        for (size_t col = 0; col < A.M; ++col) {
            A[row][col] = index(col, row);
        }
    }
    return A;
}

Matrix Matrix::copy() const {
    Matrix A(N, M);
    A._data.assign(data.begin(), data.end());
    return A;
}

bool Matrix::operator==(const Matrix& A) const {
    compatibility_check(A);
    Matrix diff = operator-(A); // difference matrix
    return std::accumulate(diff.begin(), diff.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
}

Row Matrix::row(int i) {return Row(_data, N, M, i);}
Row Matrix::operator[](int i) {return Row(_data, N, M, i);}
Column Matrix::col(int i) {return Column(_data, N, M, i);}

const ConstRow Matrix::row(int i) const {return ConstRow(data, N, M, i);}
const ConstRow Matrix::operator[](int i) const {return ConstRow(data, N, M, i);}
const ConstColumn Matrix::col(int i) const {return ConstColumn(data, N, M, i);}