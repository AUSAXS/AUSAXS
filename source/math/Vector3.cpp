#include "math/Vector3.h"
#include "math/Matrix.h"

void Vector3::rotate(Matrix& matrix) {
    *this = matrix*(*this);
}

void Vector3::rotate(Vector3& axis, const double& angle) {
    // we use the Euler-Rodrigues formulation
    axis.normalize();
    double a = cos(angle/2);
    double b = sin(angle/2);
    double c = b;
    double d = b;
    b *= axis.x;
    c *= axis.y;
    d *= axis.z;

    double aa = a*a, bb = b*b, cc = c*c, dd = d*d;
    double bc = b*c, ad = a*d, ac = a*c, ab = a*b, bd = b*d, cd = c*d;

    Matrix R{{aa+bb-cc-dd, 2*(bc-ad),   2*(bd+ac)}, 
            {2*(bc+ad),   aa+cc-bb-dd, 2*(cd-ab)},
            {2*(bd-ac),   2*(cd+ab),   aa+dd-bb-cc}};

    rotate(R);
}
