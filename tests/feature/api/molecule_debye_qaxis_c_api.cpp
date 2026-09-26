// The q values the C API returns alongside a profile must be the ones the profile was computed at, also when qmin is raised.

#include <api/api_pyausaxs.h>
#include <settings/HistogramSettings.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using namespace ausaxs;

TEST_CASE("molecule_debye: q axis honours qmin") {
    // a small deterministic cluster; the Simple path is enough since only the q labelling is under test
    std::vector<double> x, y, z, w;
    for (int i = 0; i < 40; ++i) {
        x.push_back(7*std::sin(1.3*i)); y.push_back(5*std::cos(0.7*i)); z.push_back(0.4*i - 8); w.push_back(6);
    }
    int status = 1;
    int mol = molecule_from_arrays(x.data(), y.data(), z.data(), w.data(), x.size(), &status);
    REQUIRE(status == 0);

    double qmin = GENERATE(1e-4, 0.05, 0.2);
    settings::axes::qmin = qmin;
    settings::axes::qmax = 0.5;

    // each profile is compared against its own _userq variant right away: a _raw call changes the global exv state (BL-138)
    double *q, *I, *q_raw, *I_raw, *q_ex, *I_ex;
    int n = 0, n_raw = 0, n_ex = 0;
    molecule_debye(mol, &q, &I, &n, &status);
    REQUIRE(status == 0);
    std::vector<double> I_user(n);
    molecule_debye_userq(mol, q, I_user.data(), n, &status);
    REQUIRE(status == 0);

    molecule_debye_raw(mol, &q_raw, &I_raw, &n_raw, &status);
    REQUIRE(status == 0);
    std::vector<double> I_raw_user(n_raw);
    molecule_debye_raw_userq(mol, q_raw, I_raw_user.data(), n_raw, &status);
    REQUIRE(status == 0);

    molecule_debye_exact(mol, &q_ex, &I_ex, &n_ex, &status);
    REQUIRE(status == 0);

    // all three families report the same q grid, starting at qmin
    REQUIRE(n == n_ex);
    REQUIRE(n_raw == n_ex);
    CHECK_THAT(q[0], Catch::Matchers::WithinAbs(qmin, 5e-3));
    for (int i = 0; i < n; ++i) {
        CHECK_THAT(q[i], Catch::Matchers::WithinAbs(q_ex[i], 1e-12));
        CHECK_THAT(q_raw[i], Catch::Matchers::WithinAbs(q_ex[i], 1e-12));
    }

    // evaluating a profile at its own q values reproduces it
    for (int i = 0; i < n; ++i) {
        CHECK_THAT(I_user[i], Catch::Matchers::WithinRel(I[i], 1e-6));
        CHECK_THAT(I_raw_user[i], Catch::Matchers::WithinRel(I_raw[i], 1e-6));
    }

    settings::axes::qmin = 1e-4;
}
