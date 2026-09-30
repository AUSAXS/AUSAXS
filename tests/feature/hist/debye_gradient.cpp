#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/DebyeGradient.h>
#include <settings/GeneralSettings.h>
#include <settings/MoleculeSettings.h>

#include <cmath>
#include <random>

using namespace ausaxs;
using namespace ausaxs::data;

// L = sum_q v(q) I(q) for the exact raw Debye sum I(q) = sum_ij w_i w_j sinc(q r_ij)
static double adjoint_loss(const std::vector<AtomFF>& atoms, const std::vector<double>& q, const std::vector<double>& v) {
    double L = 0;
    for (int k = 0; k < static_cast<int>(q.size()); ++k) {
        double I = 0;
        for (const auto& a : atoms) {
            for (const auto& b : atoms) {
                double qr = q[k]*a.coordinates().distance(b.coordinates());
                I += a.weight()*b.weight()*(qr < 1e-12 ? 1 : std::sin(qr)/qr);
            }
        }
        L += v[k]*I;
    }
    return L;
}

// the analytic gradient evaluated directly, without the tabulation of the pair derivative
static std::vector<Vector3<double>> direct_gradient(const std::vector<AtomFF>& atoms, const std::vector<double>& q, const std::vector<double>& v) {
    std::vector<Vector3<double>> g(atoms.size(), {0, 0, 0});
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        for (int j = 0; j < static_cast<int>(atoms.size()); ++j) {
            if (i == j) {continue;}
            Vector3<double> d = atoms[i].coordinates() - atoms[j].coordinates();
            double r = d.norm();
            double dsinc = 0;
            for (int k = 0; k < static_cast<int>(q.size()); ++k) {
                double qr = q[k]*r;
                dsinc += v[k]*(std::cos(qr) - std::sin(qr)/qr)/r;
            }
            g[i] += d*(2*atoms[i].weight()*atoms[j].weight()*dsinc/r);
        }
    }
    return g;
}

static double max_relative_deviation(const std::vector<Vector3<double>>& a, const std::vector<Vector3<double>>& b) {
    double scale = 0, deviation = 0;
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        scale = std::max(scale, b[i].norm());
        deviation = std::max(deviation, (a[i] - b[i]).norm());
    }
    return deviation/scale;
}

TEST_CASE("debye_raw_vjp: agrees with finite differences") {
    settings::general::verbose = false;
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> coordinate(-15, 15);
    std::uniform_real_distribution<double> weight(1, 8), adjoint(-1, 1);

    std::vector<AtomFF> atoms;
    for (int i = 0; i < 40; ++i) {
        atoms.emplace_back(Vector3<double>(coordinate(gen), coordinate(gen), coordinate(gen)), form_factor::form_factor_t::C);
        atoms.back().weight() = weight(gen);
    }
    std::vector<double> q, v;
    for (int k = 0; k < 25; ++k) {q.push_back(0.02*k); v.push_back(adjoint(gen));}

    Molecule molecule({Body{atoms}});
    auto g = hist::debye_raw_vjp(molecule, q, v);

    constexpr double h = 1e-4;
    std::vector<Vector3<double>> fd(atoms.size());
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        for (int c = 0; c < 3; ++c) {
            auto plus = atoms, minus = atoms;
            plus[i].coordinates()[c] += h;
            minus[i].coordinates()[c] -= h;
            fd[i][c] = (adjoint_loss(plus, q, v) - adjoint_loss(minus, q, v))/(2*h);
        }
    }
    REQUIRE(max_relative_deviation(g, fd) < 1e-5);
}

TEST_CASE("debye_raw_vjp: agrees with the direct analytic gradient") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    Molecule molecule("tests/files/2epe.pdb");
    molecule.clear_hydration();

    std::mt19937 gen(7);
    std::uniform_real_distribution<double> adjoint(-1, 1);
    std::vector<double> q, v;
    for (int k = 0; k < 20; ++k) {q.push_back(0.001 + 0.05*k); v.push_back(adjoint(gen));}

    auto g = hist::debye_raw_vjp(molecule, q, v);
    auto direct = direct_gradient(molecule.get_atoms(), q, v);
    REQUIRE(max_relative_deviation(g, direct) < 1e-5);
}

TEST_CASE("debye_raw_vjp: rejects hydrated molecules") {
    settings::general::verbose = false;
    Molecule molecule("tests/files/2epe.pdb");
    molecule.generate_new_hydration();
    REQUIRE_THROWS(hist::debye_raw_vjp(molecule, {0.1}, {1}));
}
