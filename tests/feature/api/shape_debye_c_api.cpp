#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <api/api_pyausaxs.h>
#include <settings/HistogramSettings.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

// the Gaussian form factor of one lattice cell of volume V = spacing^3, as in the Grid excluded volume model
static double cell_ff(double q, double spacing) {
    double V = std::pow(spacing, 3);
    return V*std::exp(-std::pow(V, 2./3)*q*q/(4*std::numbers::pi));
}

TEST_CASE("shape_debye_userq: two cells") {
    // one pair distance, so the weighted bin centre is the distance itself and the result is exact
    for (double spacing : {1.0, 0.5}) {
        std::vector<double> x = {0, 3*spacing}, y = {spacing, spacing}, z = {-2*spacing, 2*spacing};
        double d = 5*spacing;
        std::vector<double> q = {0, 0.1, 0.3, 0.7}, I(q.size());
        int status = 1;
        shape_debye_userq(x.data(), y.data(), z.data(), 2, spacing, q.data(), I.data(), q.size(), &status);
        REQUIRE(status == 0);
        for (unsigned int i = 0; i < q.size(); ++i) {
            double sinc = q[i] == 0 ? 1 : std::sin(q[i]*d)/(q[i]*d);
            CHECK_THAT(I[i], Catch::Matchers::WithinRel((2 + 2*sinc)*std::pow(cell_ff(q[i], spacing), 2), 1e-5));
        }
    }
}

TEST_CASE("shape_debye_userq: homogeneous sphere") {
    // the grid models a homogeneous sphere, so it must reproduce the analytical intensity of a sphere of the same volume up
    // to the boundary effects of the lattice representation (~4% of the local curve height by q = 0.5); measured against
    // the local curve height, since the relative error is meaningless at its zeros
    ausaxs::settings::axes::bin_width = 0.05; // the default bins are the dominant error for a body this smooth
    double R = 20, spacing = 1;
    std::vector<double> x, y, z;
    for (int i = -21; i <= 21; ++i) {
        for (int j = -21; j <= 21; ++j) {
            for (int k = -21; k <= 21; ++k) {
                if (R*R < i*i + j*j + k*k) {continue;}
                x.push_back(i + 7); y.push_back(j - 3); z.push_back(k); // an arbitrary offset
            }
        }
    }
    int n = x.size();
    double R_eff = std::cbrt(3.*n/(4*std::numbers::pi));

    std::vector<double> q, I;
    for (double qi = 1e-3; qi <= 0.5; qi += 0.005) {q.push_back(qi);}
    I.resize(q.size());
    int status = 1;
    shape_debye_userq(x.data(), y.data(), z.data(), n, spacing, q.data(), I.data(), q.size(), &status);
    REQUIRE(status == 0);
    CHECK_THAT(I[0], Catch::Matchers::WithinRel(std::pow(n*cell_ff(q[0], spacing), 2), 1e-4));

    std::vector<double> exact(q.size());
    for (unsigned int i = 0; i < q.size(); ++i) {
        double s = q[i]*R_eff, V = n*std::pow(spacing, 3);
        exact[i] = std::pow(V*3*(std::sin(s) - s*std::cos(s))/(s*s*s), 2);
    }
    double height = 0;
    for (int i = static_cast<int>(q.size())-1; 0 <= i; --i) {
        height = std::max(height, exact[i]);
        CHECK(std::abs(I[i] - exact[i]) < 5e-2*height);
    }
    ausaxs::settings::axes::bin_width = 0.25;
}

TEST_CASE("shape_debye_userq: invalid input") {
    std::vector<double> q = {0.1}, I(1);
    int status = 0;
    std::vector<double> x = {0, 0}, y = {0, 0}, z = {0, 1.5};
    shape_debye_userq(x.data(), y.data(), z.data(), 2, 1, q.data(), I.data(), 1, &status);
    CHECK(status != 0); // off the lattice

    status = 0;
    z = {1, 1};
    shape_debye_userq(x.data(), y.data(), z.data(), 2, 1, q.data(), I.data(), 1, &status);
    CHECK(status != 0); // the same site twice

    status = 0;
    shape_debye_userq(x.data(), y.data(), z.data(), 0, 1, q.data(), I.data(), 1, &status);
    CHECK(status != 0); // empty
}
