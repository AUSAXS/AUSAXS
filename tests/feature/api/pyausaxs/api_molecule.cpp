// Tests of the molecule_* functions of the C API backing pyausaxs.

#include <api/ObjectStorage.h>
#include <api/api_pyausaxs.h>
#include <api/pyausaxs/api_settings.h>
#include <constants/ConstantsAxes.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/atoms/AtomFF.h>
#include <form_factor/xray/XrayFormFactor.h>
#include <math/ConstexprMath.h>
#include <settings/HistogramSettings.h>
#include <utility/Random.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <random>
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

    // each profile is compared against its own _userq variant
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
}

TEST_CASE("molecule_distance_histogram_raw: weight products without excluded volume") {
    // signed weights, so a subtracted excluded volume or a dropped self-term cannot cancel out
    std::vector<double> x, y, z, w;
    double w_sum = 0;
    for (int i = 0; i < 60; ++i) {
        x.emplace_back(9*std::sin(1.7*i)); y.emplace_back(6*std::cos(0.9*i)); z.emplace_back(0.5*i - 15); w.emplace_back(1.5*std::sin(2.3*i) + 0.5);
        w_sum += w.back();
    }
    int status = 1;
    int mol = molecule_from_arrays(x.data(), y.data(), z.data(), w.data(), static_cast<int>(x.size()), &status);
    REQUIRE(status == 0);

    // Fraser cannot evaluate these atoms at all since they carry no form factor, but the raw histogram never uses it
    std::string exv = GENERATE("simple", "none", "fraser");
    set_setting("exv_model", exv.c_str(), &status);
    REQUIRE(status == 0);

    double *counts, *axis;
    int n_bins = 0;
    molecule_distance_histogram_raw(mol, &counts, &axis, &n_bins, &status);
    REQUIRE(status == 0);
    REQUIRE(n_bins > 0);

    double total = 0;
    for (int i = 0; i < n_bins; ++i) {total += counts[i];}
    CHECK_THAT(total, Catch::Matchers::WithinRel(w_sum*w_sum, 1e-6)); // the histogram kernels take the weights in single precision
    CHECK(axis[0] == 0);
}

// raw profiles must hold at any q: a protein with random weights, compared with the direct sum up to q = 10. Proteins are disordered
// enough for unweighted bins in the SAXS range, but not this far beyond it.
TEST_CASE("molecule_debye_raw_userq: accurate at high q") {
    int status = 1;
    int pdb = pdb_read("tests/files/2epe.pdb", &status);
    REQUIRE(status == 0);
    int *serial, *res_seq;
    const char **name, **alt_name, **res_name, **i_code, **element, **charge;
    const char* chain_id;
    double *px, *py, *pz, *occupancy, *b_factor;
    int n_atoms = 0;
    pdb_get_data(pdb, &serial, &name, &alt_name, &res_name, &chain_id, &res_seq, &i_code, &px, &py, &pz, &occupancy, &b_factor, &element, &charge, &n_atoms, &status);
    REQUIRE(status == 0);

    auto& rng = random::generator();
    std::uniform_real_distribution<double> weight(0, 1);
    std::vector<double> x(px, px + n_atoms), y(py, py + n_atoms), z(pz, pz + n_atoms), w(n_atoms);
    for (auto& wi : w) {wi = weight(rng);}
    int mol = molecule_from_arrays(x.data(), y.data(), z.data(), w.data(), n_atoms, &status);
    REQUIRE(status == 0);
    settings::hist::weighted_bins = settings::hist::WeightedBins::Value::Auto;

    std::vector<double> q(100), I(q.size());
    for (unsigned int k = 0; k < q.size(); ++k) {q[k] = 0.1 + 0.1*k;}
    molecule_debye_raw_userq(mol, q.data(), I.data(), static_cast<int>(q.size()), &status);
    REQUIRE(status == 0);
    CHECK(settings::hist::weighted_bins->is_auto());

    std::vector<double> expected(q.size());
    molecule_debye_exact_userq(mol, q.data(), expected.data(), static_cast<int>(q.size()), &status);
    REQUIRE(status == 0);
    for (unsigned int k = 0; k < q.size(); ++k) {
        CHECK_THAT(I[k], Catch::Matchers::WithinRel(expected[k], 0.03)); // the tolerance of SasView's own comparison
    }
}

// atoms on a line at the marks of a Golomb ruler, so every distance is distinct and has a weighted bin of its own, and the binned sum must
// match the direct one
TEST_CASE("molecule_debye_raw_ff: agrees with the direct sum with atomic form factors") {
    using ff_t = form_factor::form_factor_t;
    constexpr std::array marks = {0, 1, 5, 12, 25, 27, 35};
    constexpr std::array types = {ff_t::C, ff_t::N, ff_t::O, ff_t::S, ff_t::Fe, ff_t::Zn, ff_t::I};
    std::vector<data::AtomFF> atoms;
    for (unsigned int i = 0; i < marks.size(); ++i) {atoms.emplace_back(Vector3<double>{1.13*marks[i], 0.2*i, 0}, types[i]);}
    int mol = api::ObjectStorage::register_object(data::Molecule({data::Body{atoms}}));

    int status = 1;
    double *q, *I;
    int n = 0;
    molecule_debye_raw_ff(mol, &q, &I, &n, &status);
    REQUIRE(status == 0);
    REQUIRE(n > 0);

    auto f = [] (ff_t t, double q) {return form_factor::xray::raw::get(t).evaluate(q);};
    for (int k = 0; k < n; ++k) {
        double expected = 0, scale = 0;
        for (unsigned int i = 0; i < atoms.size(); ++i) {
            expected += form_factor::xray::raw::get(types[i]).evaluate_self(q[k]);
            scale += std::abs(f(types[i], q[k]));
            for (unsigned int j = 0; j < atoms.size(); ++j) {
                if (i == j) {continue;}
                expected += f(types[i], q[k])*f(types[j], q[k])*math::sinc(q[k]*atoms[i].get_atom().coordinates().distance(atoms[j].get_atom().coordinates()));
            }
        }
        CHECK_THAT(I[k], Catch::Matchers::WithinAbs(expected, 1e-5*scale*scale)); // the intensity crosses zero, so measure against its scale
    }

    SECTION("userq reproduces the profile, but refuses q beyond the default axis") {
        std::vector<double> qv(q, q + n), Iv(n);
        molecule_debye_raw_ff_userq(mol, qv.data(), Iv.data(), n, &status);
        REQUIRE(status == 0);
        for (int k = 0; k < n; ++k) {CHECK_THAT(Iv[k], Catch::Matchers::WithinRel(I[k], 1e-9));}

        std::array q_high = {0.5, constants::axes::q_axis.max + 1};
        std::array<double, 2> I_high{};
        molecule_debye_raw_ff_userq(mol, q_high.data(), I_high.data(), 2, &status);
        CHECK(status != 0);
    }
}

TEST_CASE("molecule_debye_raw_ff: molecules without form factors are refused") {
    std::array x = {0.0, 1.0, 0.0};
    std::array y = {0.0, 0.0, 1.0};
    std::array z = {0.0, 0.0, 0.0};
    std::array w = {1.0, 1.0, 1.0};
    int status = 1;
    int mol = molecule_from_arrays(x.data(), y.data(), z.data(), w.data(), 3, &status);
    REQUIRE(status == 0);

    double *q, *I;
    int n = 0;
    molecule_debye_raw_ff(mol, &q, &I, &n, &status);
    CHECK(status != 0);
}

TEST_CASE("object ids: an object of the wrong type is refused") {
    // an id is just an int from the caller, so it can name any stored object; reading it as the wrong type aborts or segfaults the host
    int status = 0;
    int mol = molecule_from_file("tests/files/2epe.pdb", &status);
    REQUIRE(status == 0);
    int pdb = pdb_read("tests/files/2epe.pdb", &status);
    REQUIRE(status == 0);
    int data = data_read("tests/files/2epe.dat", &status);
    REQUIRE(status == 0);

    int fit = molecule_debye_fit(mol, data, &status);
    REQUIRE(status == 0); // control

    for (int not_data : {mol, pdb, fit}) {
        molecule_debye_fit(mol, not_data, &status);
        CHECK(status != 0);
    }
    for (int not_molecule : {data, pdb, fit}) {
        molecule_debye_fit(not_molecule, data, &status);
        CHECK(status != 0);
        iterative_fit_init(not_molecule, &status);
        CHECK(status != 0);
        double *q = nullptr, *I = nullptr;
        int n = 0;
        molecule_debye_exact(not_molecule, &q, &I, &n, &status);
        CHECK(status != 0);
    }
    molecule_from_pdb_id(mol, &status);
    CHECK(status != 0);
}

TEST_CASE("object ids: a missing id and a wrong type give different errors") {
    auto last_error = [] {
        char* msg = nullptr;
        int status = 0;
        get_last_error_msg(&msg, &status);
        return std::string(msg);
    };

    int status = 0;
    int data = data_read("tests/files/2epe.dat", &status);
    REQUIRE(status == 0);
    int mol = molecule_from_file("tests/files/2epe.pdb", &status);
    REQUIRE(status == 0);

    molecule_debye_fit(mol, mol, &status);
    REQUIRE(status != 0);
    CHECK(last_error() == "Object \"" + std::to_string(mol) + "\" has the wrong type (expected: dataset).");

    deallocate(data, &status);
    REQUIRE(status == 0);
    molecule_debye_fit(mol, data, &status);
    REQUIRE(status != 0);
    CHECK(last_error() == "No object with id \"" + std::to_string(data) + "\" exists (expected: dataset). It may have been deallocated.");
}
