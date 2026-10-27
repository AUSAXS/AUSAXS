#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <api/api_sasview.h>
#include <settings/ExvSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/MoleculeSettings.h>
#include <utility/Random.h>

#include <cmath>
#include <vector>

using namespace ausaxs;

// the plain weighted Debye sum, evaluated directly
static std::vector<double> debye_exact(const std::vector<double>& q, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& w) {
    std::vector<double> I(q.size(), 0);
    for (unsigned int k = 0; k < q.size(); ++k) {
        for (unsigned int i = 0; i < x.size(); ++i) {
            for (unsigned int j = 0; j < x.size(); ++j) {
                double r = std::sqrt(std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2) + std::pow(z[i]-z[j], 2));
                double qr = q[k]*r;
                I[k] += w[i]*w[j]*(qr < 1e-12 ? 1 : std::sin(qr)/qr);
            }
        }
    }
    return I;
}

TEST_CASE("debye_no_ff: matches the direct Debye sum") {
    auto& rng = random::generator();
    std::uniform_real_distribution<double> pos(-20, 20), weight(0.5, 2);
    int N = 50;
    std::vector<double> x(N), y(N), z(N), w(N);
    for (int i = 0; i < N; ++i) {x[i] = pos(rng); y[i] = pos(rng); z[i] = pos(rng); w[i] = weight(rng);}
    std::vector<double> q = {0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5}, I(q.size());

    int status = -1;
    debye_no_ff(q.data(), x.data(), y.data(), z.data(), w.data(), static_cast<int>(q.size()), N, I.data(), &status);
    REQUIRE(status == 0);
    auto expected = debye_exact(q, x, y, z, w);
    for (unsigned int k = 0; k < q.size(); ++k) {
        // measured against I(0), since the binned distances shift the curve slightly where it is small (worst seen: 3e-4 over 300 random draws)
        CHECK_THAT(I[k], Catch::Matchers::WithinAbs(expected[k], 1e-3*expected[0]));
    }
}

// a SasView call must not change the settings of later, unrelated calculations in the same process
TEST_CASE("debye_no_ff: leaves the caller's settings untouched") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
    settings::molecule::implicit_hydrogens = true;
    settings::axes::qmax = 0.4;

    std::vector<double> x = {0, 0}, y = {0, 0}, z = {0, 5}, w = {1, 1}, q = {0, 0.1, 0.3}, I(q.size());
    int status = -1;
    debye_no_ff(q.data(), x.data(), y.data(), z.data(), w.data(), static_cast<int>(q.size()), 2, I.data(), &status);
    REQUIRE(status == 0);
    CHECK_THAT(I[1], Catch::Matchers::WithinRel(2 + 2*std::sin(0.5)/0.5, 1e-6));

    CHECK(settings::exv::exv_method == settings::exv::ExvMethod::Fraser);
    CHECK(settings::molecule::implicit_hydrogens == true);
    CHECK(settings::axes::qmax == 0.4);
}
