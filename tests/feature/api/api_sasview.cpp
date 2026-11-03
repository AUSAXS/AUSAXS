#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <api/api_sasview.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/scattering_model/ExactDebyeCalculator.h>
#include <settings/ExvSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/MoleculeSettings.h>
#include <utility/Random.h>

#include <cmath>
#include <random>
#include <vector>

using namespace ausaxs;

// the built-in exact Debye sum over the same weighted points, without the unit form factor it applies
static std::vector<double> debye_exact(const std::vector<double>& q, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& w) {
    std::vector<data::Atom> atoms;
    atoms.reserve(x.size());
    for (unsigned int i = 0; i < x.size(); ++i) {atoms.emplace_back(Vector3<double>{x[i], y[i], z[i]}, w[i]);}
    auto I = hist::exact_debye_transform(data::Molecule({data::Body{atoms}}), q);
    for (unsigned int k = 0; k < q.size(); ++k) {I[k] *= std::exp(q[k]*q[k]);}
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

// SasView's own check of this function: a protein with random weights, compared with the direct sum up to q = 10. Proteins are disordered
// enough for unweighted bins in the SAXS range, but not this far beyond it.
TEST_CASE("debye_no_ff: accurate at high q") {
    settings::hist::weighted_bins = settings::hist::WeightedBins::Value::Auto;

    data::Molecule protein("tests/files/2epe.pdb");
    auto& rng = random::generator();
    std::uniform_real_distribution<double> weight(0, 1);
    std::vector<double> x, y, z, w;
    for (const auto& atom : protein.iterate_atoms()) {
        x.emplace_back(atom.coordinates().x()); y.emplace_back(atom.coordinates().y()); z.emplace_back(atom.coordinates().z());
        w.emplace_back(weight(rng));
    }
    std::vector<double> q(100), I(q.size());
    for (unsigned int k = 0; k < q.size(); ++k) {q[k] = 0.1 + 0.1*k;}

    int status = -1;
    debye_no_ff(q.data(), x.data(), y.data(), z.data(), w.data(), static_cast<int>(q.size()), static_cast<int>(x.size()), I.data(), &status);
    REQUIRE(status == 0);
    auto expected = debye_exact(q, x, y, z, w);
    for (unsigned int k = 0; k < q.size(); ++k) {
        CHECK_THAT(I[k], Catch::Matchers::WithinRel(expected[k], 0.03)); // the tolerance of SasView's own comparison
    }
}

// a SasView call must not change the settings of later, unrelated calculations in the same process
TEST_CASE("debye_no_ff: leaves the caller's settings untouched") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
    settings::molecule::implicit_hydrogens = true;
    settings::axes::qmax = 0.4;
    settings::hist::weighted_bins = settings::hist::WeightedBins::Value::Auto;

    std::vector<double> x = {0, 0}, y = {0, 0}, z = {0, 5}, w = {1, 1}, q = {0, 0.1, 0.3}, I(q.size());
    int status = -1;
    debye_no_ff(q.data(), x.data(), y.data(), z.data(), w.data(), static_cast<int>(q.size()), 2, I.data(), &status);
    REQUIRE(status == 0);
    CHECK_THAT(I[1], Catch::Matchers::WithinRel(2 + 2*std::sin(0.5)/0.5, 1e-6));

    CHECK(settings::exv::exv_method == settings::exv::ExvMethod::Fraser);
    CHECK(settings::molecule::implicit_hydrogens == true);
    CHECK(settings::axes::qmax == 0.4);
    CHECK(settings::hist::weighted_bins->is_auto());
}
