#include <crystal/Fval.h>

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

std::vector<Vector3<double>>& Fval::get_points() {
    return points;
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
