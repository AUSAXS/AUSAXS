#include <crystal/Fval.h>
#include <Symbols.h>
#include <data/Protein.h>

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

Protein Fval::as_protein() {
    std::vector<Atom> atoms(points.size());
    for (unsigned int i = 0; i < points.size(); i++) {
        atoms[i] = Atom(points[i], 1, "H", "LYS", i);
    } 
    return Protein(atoms);
}

// #include <settings/CrystalSettings.h>
// std::vector<std::complex<double>> x_factors;
// std::vector<std::complex<double>> y_factors;
// std::vector<std::complex<double>> z_factors;
// std::vector<std::complex<double>> int_factors;
// void Fval::precompute_factors() {
//     x_factors.reserve(points.size());
//     y_factors.reserve(points.size());
//     z_factors.reserve(points.size());

//     for (const Vector3<double>& point : points) {
//         std::complex<double> x_factor = std::polar(1.0, -2*M_PI*ap.x()*point.x());
//         std::complex<double> y_factor = std::polar(1.0, -2*M_PI*bp.y()*point.y());
//         std::complex<double> z_factor = std::polar(1.0, -2*M_PI*cp.z()*point.z());

//         x_factors.push_back(x_factor);
//         y_factors.push_back(y_factor);
//         z_factors.push_back(z_factor);
//     }

//     auto max = std::max({settings::crystal::h, settings::crystal::k, settings::crystal::l});
//     int_factors.reserve(2*max+1);
//     for (double i = -max; i <= max; ++i) {
//         int_factors.push_back(std::polar(1.0, i));
//     }
// }

std::complex<double> Fval::F() {
    std::complex<double> result = 0;
    for (const Vector3<double>& point : points) {
        result += std::polar(1.0, -q.dot(point));
    }
    return result;
}

Vector3<double> Fval::Q() const {
    return Q(hkl);
}

Vector3<double> Fval::Q(const Miller& miller) {
    return miller.h*ap + miller.k*bp + miller.l*cp;
}