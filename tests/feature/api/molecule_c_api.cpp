// Tests of the molecule_* functions of the C API backing pyausaxs.

#include <api/api_pyausaxs.h>
#include <api/pyausaxs/api_settings.h>
#include <settings/HistogramSettings.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace ausaxs;

TEST_CASE("molecule_distance_histogram: UNKNOWN form factors with Fraser exv model") {
    // Create atoms without form factor information (like molecule_from_arrays does)
    std::array x = {0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, -1.0};
    std::array y = {0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0, -1.0};
    std::array z = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 1.0, -1.0};
    std::array w = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int n_atoms = 9;

    int status = 0;

    // Create molecule from arrays
    int mol_id = molecule_from_arrays(x.data(), y.data(), z.data(), w.data(), n_atoms, &status);
    REQUIRE(status == 0);
    REQUIRE(mol_id >= 0);

    // Set problematic settings via C API
    set_setting("exv_model", "Fraser", &status);
    REQUIRE(status == 0);

    set_setting("N", "10", &status);
    set_setting("excluded_volume", "true", &status);
    set_setting("solvent_density", "true", &status);
    REQUIRE(status == 0);

    // This should fail when trying to create the histogram because Fraser ExV requires form factor info
    double *aa, *aw, *ww, *axis;
    int n_bins;
    molecule_distance_histogram(mol_id, &aa, &aw, &ww, &axis, &n_bins, &status);

    // We expect this to fail because the Fraser model requires form factor information
    REQUIRE(status != 0);

    char* error_msg = nullptr;
    int error_status = 0;
    get_last_error_msg(&error_msg, &error_status);
    REQUIRE(error_msg != nullptr);

    std::string error_str(error_msg);
    CHECK(error_str.find("UNKNOWN form factor") != std::string::npos);
}

// The q values returned alongside a profile must be the ones the profile was computed at, also when qmin is raised.
TEST_CASE("molecule_debye: q axis honours qmin") {
    // a small deterministic cluster; the Simple path is enough since only the q labelling is under test
    std::vector<double> x, y, z, w;
    for (int i = 0; i < 40; ++i) {
        x.push_back(7*std::sin(1.3*i)); y.push_back(5*std::cos(0.7*i)); z.push_back(0.4*i - 8); w.push_back(6);
    }
    int status = 1;
    int mol = molecule_from_arrays(x.data(), y.data(), z.data(), w.data(), static_cast<int>(x.size()), &status);
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
