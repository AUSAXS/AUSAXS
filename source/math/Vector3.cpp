#include "math/Vector3.h"
#include "math/Matrix.h"

void Vector3::rotate(const Matrix& matrix) {
    *this = matrix*(*this);
}

void Vector3::rotate(const Vector3& axis, const double angle) {
    // we use the Euler-Rodrigues formulation
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

    rotate(R);
}
