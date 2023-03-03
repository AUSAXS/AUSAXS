#include <crystal/Fval.h>
#include <Symbols.h>

#include <cmath>

using namespace crystal;

Fval::Fval() = default;
Fval::Fval(double h, double k, double l) : hkl(h, k, l) {
    q = Q();
    qlength = q.norm();
    fval = F();
}

void Fval::set_points(std::vector<Vector3<double>>&& points) {
    Fval::points = std::move(points);
}

void Fval::set_basis(const Basis3D& basis) {
    ap = {2*M_PI/basis.x.x(), 0, 0};
    bp = {0, 2*M_PI/basis.y.y(), 0};
    cp = {0, 0, 2*M_PI/basis.z.z()};
}

std::vector<Vector3<double>>& Fval::get_points() {
    return points;
}

Basis3D Fval::get_basis() {
    return Basis3D(ap, bp, cp);
}

double Fval::I() const {
    return std::norm(fval);
}

std::complex<double> Fval::F() {
    std::complex<double> result = 0;
    for (const Vector3<double>& point : points) {
        result += exp(-i*q.dot(point));
    }
    return result;
}

Vector3<double> Fval::Q() const {
    return hkl.h*ap + hkl.k*bp + hkl.l*cp;
}
