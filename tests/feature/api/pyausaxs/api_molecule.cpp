// Tests of the molecule_* functions of the C API backing pyausaxs.

#include <api/ObjectStorage.h>
#include <api/api_pyausaxs.h>
#include <api/pyausaxs/api_settings.h>
#include <data/Molecule.h>
#include <hist/ScatteringProfile.h>
#include <hist/scattering_model/IScatteringModel.h>
#include <hist/scattering_model/amplitude/ScatteringAmplitude.h>
#include <math/Vector3.h>
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

// end-to-end: with non-default settings, the profile functions without q values lie on get_q_axis, and each returns the profile of the
// molecule's own scattering model or amplitude. The values themselves are verified by the tests of the core library.
TEST_CASE("profile functions with user settings") {
    int status = 1;
    int mol = molecule_from_file("tests/files/2epe.pdb", &status);
    REQUIRE(status == 0);
    molecule_hydrate(mol, &status);
    REQUIRE(status == 0);
    for (auto [key, value] : std::array<std::array<const char*, 2>, 3>{{{"exv_model", "fraser"}, {"qmin", "0.05"}, {"qmax", "0.4"}}}) {
        set_setting(key, value, &status);
        REQUIRE(status == 0);
    }

    // the references are taken first, since the raw profiles rebuild the histogram manager with their own settings
    std::array<double, 3> centre = {50, 0, 0};
    auto& molecule = api::ObjectStorage::get<data::Molecule>(mol, "molecule");
    auto I_ref = molecule.get_histogram()->debye_transform();
    auto q_ref = I_ref.get_axis().as_vector();
    auto A_centroid_ref = molecule.get_scattering_amplitude()->evaluate(q_ref);
    auto A_centre_ref = molecule.get_scattering_amplitude()->evaluate(q_ref, Vector3<double>(centre[0], centre[1], centre[2]));

    const double* q;
    int n = 0;
    get_q_axis(&q, &n, &status);
    REQUIRE(status == 0);
    REQUIRE(n == static_cast<int>(q_ref.size()));
    for (int k = 0; k < n; ++k) {CHECK_THAT(q[k], Catch::Matchers::WithinRel(q_ref[k], 1e-12));}
    std::vector<double> q_user(q, q + n);

    auto check = [n] (const double* profile, int n_profile, const std::vector<double>& ref) {
        REQUIRE(n_profile == n);
        for (int k = 0; k < n; ++k) {CHECK_THAT(profile[k], Catch::Matchers::WithinRel(ref[k], 1e-6));}
    };

    double *I, *A, *I_beta, *A_beta;
    int n_out = 0;
    molecule_debye(mol, &I, &n_out, &status);
    REQUIRE(status == 0);
    check(I, n_out, I_ref.get_intensity());

    molecule_amplitude_centroid(mol, &A, &n_out, &status);
    REQUIRE(status == 0);
    check(A, n_out, A_centroid_ref);
    std::vector<double> A_user(n);
    molecule_amplitude_centroid_userq(mol, q_user.data(), A_user.data(), n, &status);
    REQUIRE(status == 0);
    check(A_user.data(), n, A_centroid_ref);

    molecule_amplitude(mol, centre.data(), &A, &n_out, &status);
    REQUIRE(status == 0);
    check(A, n_out, A_centre_ref);
    molecule_amplitude_userq(mol, centre.data(), q_user.data(), A_user.data(), n, &status);
    REQUIRE(status == 0);
    check(A_user.data(), n, A_centre_ref);

    molecule_beta_centroid(mol, &I_beta, &A_beta, &n_out, &status);
    REQUIRE(status == 0);
    check(I_beta, n_out, I_ref.get_intensity());
    check(A_beta, n_out, A_centroid_ref);
    molecule_beta(mol, centre.data(), &I_beta, &A_beta, &n_out, &status);
    REQUIRE(status == 0);
    check(I_beta, n_out, I_ref.get_intensity());
    check(A_beta, n_out, A_centre_ref);

    // the raw and exact profiles follow other paths; only their axis is checked here
    molecule_debye_raw(mol, &I, &n_out, &status);
    REQUIRE(status == 0);
    CHECK(n_out == n);
    molecule_debye_exact(mol, &I, &n_out, &status);
    REQUIRE(status == 0);
    CHECK(n_out == n);
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
        double* I = nullptr;
        int n = 0;
        molecule_debye_exact(not_molecule, &I, &n, &status);
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
