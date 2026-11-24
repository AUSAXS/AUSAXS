#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/scattering_model/IScatteringModel.h>
#include <hist/scattering_model/IScatteringModelExv.h>
#include <hist/scattering_model/amplitude/ScatteringAmplitude.h>
#include <math/ConstexprMath.h>
#include <settings/All.h>

#include <cmath>

using namespace ausaxs;
using namespace ausaxs::data;

namespace {
    // <A>^2 <= <A^2> at every q, with equality at q = 0, where every phase is 1
    void check_consistent(const std::vector<double>& q, const std::vector<double>& A, const std::vector<double>& I) {
        REQUIRE(A.size() == I.size());
        CHECK_THAT(A.front()*A.front()/I.front(), Catch::Matchers::WithinRel(1, 1e-3));
        for (std::size_t i = 0; i < I.size(); ++i) {
            INFO("q = " << q[i] << ", beta = " << A[i]*A[i]/I[i]);
            REQUIRE(A[i]*A[i] <= 1.01*I[i]);
        }
    }
}

// two point scatterers of weights 1 and 3 a distance d apart: the centre is 3d/4 from the first, and each bin holds a single distance, so the
// binned sums are exact
TEST_CASE("Molecule::get_scattering_amplitude: agrees with the analytical result for two points") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::None;
    constexpr double d = 10;
    Molecule molecule({Body(std::vector<AtomFF>{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C, 1),
        AtomFF({d, 0, 0}, form_factor::form_factor_t::C, 3)
    })});

    std::vector<double> q = {0, 0.01, 0.1, 0.5, 1, 2, 5, 10};
    auto A = molecule.get_scattering_amplitude()->evaluate(q);
    for (std::size_t i = 0; i < q.size(); ++i) {
        double expected = std::exp(-0.5*q[i]*q[i])*(math::sinc(q[i]*0.75*d) + 3*math::sinc(q[i]*0.25*d));
        REQUIRE_THAT(A[i], Catch::Matchers::WithinAbs(expected, 1e-6));
    }
}

TEST_CASE("Molecule::get_scattering_amplitude: the binned sums agree with the exact ones") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::None;
    Molecule molecule("tests/files/2epe.pdb");
    molecule.generate_new_hydration();
    auto centre = molecule.get_cm();

    std::vector<double> q = {0.01, 0.1, 0.2, 0.5, 1, 2, 5, 10};
    auto A = molecule.get_scattering_amplitude()->evaluate(q, centre);
    for (std::size_t i = 0; i < q.size(); ++i) {
        double exact = 0;
        for (const auto& a : molecule.iterate_atoms()) {exact += a.weight()*math::sinc(q[i]*a.coordinates().distance(centre));}
        for (const auto& w : molecule.iterate_waters()) {exact += w.weight()*math::sinc(q[i]*w.coordinates().distance(centre));}
        exact *= std::exp(-0.5*q[i]*q[i]);
        INFO("q = " << q[i]);
        REQUIRE_THAT(A[i], Catch::Matchers::WithinAbs(exact, 1e-4*A[0]));
    }
}

// the amplitude is the counterpart of the intensity of the same model
TEST_CASE("Molecule::get_scattering_amplitude: consistent with the intensity") {
    settings::general::verbose = false;
    auto exv = GENERATE(
        settings::exv::ExvMethod::Simple, settings::exv::ExvMethod::None,
        settings::exv::ExvMethod::Vacuum, settings::exv::ExvMethod::Fraser
    );
    settings::exv::exv_method = exv;
    auto partial = GENERATE(false, true);
    SECTION("exv model " + std::to_string(static_cast<int>(exv)) + (partial ? ", partial manager" : "")) {
        Molecule molecule("tests/files/2epe.pdb");
        molecule.generate_new_hydration();
        if (partial) {
            molecule.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
        }

        auto model = molecule.get_histogram();
        auto I = model->debye_transform();
        auto q = I.get_axis().as_vector();
        auto A = molecule.get_scattering_amplitude()->evaluate(q);
        check_consistent(q, A, I.get_intensity());

        // the parameters of a freshly built model are the defaults of ScatteringAmplitude::evaluate
        auto A_model = molecule.get_scattering_amplitude()->evaluate(q, std::nullopt, model->get_parameters());
        for (std::size_t i = 0; i < q.size(); ++i) {REQUIRE_THAT(A_model[i], Catch::Matchers::WithinRel(A[i], 1e-12));}
    }
}

// the amplitude uses the parameters of the model it is evaluated with
TEST_CASE("ScatteringAmplitude::evaluate: follows the fitted parameters") {
    settings::general::verbose = false;
    Molecule molecule("tests/files/2epe.pdb");
    molecule.generate_new_hydration();

    SECTION("simple") {
        settings::exv::exv_method = settings::exv::ExvMethod::Simple;
        auto model = molecule.get_histogram();
        model->apply_water_scaling_factor(3);
        auto I = model->debye_transform();
        auto q = I.get_axis().as_vector();
        auto A = molecule.get_scattering_amplitude()->evaluate(q, std::nullopt, model->get_parameters());
        check_consistent(q, A, I.get_intensity());
    }

    SECTION("fraser") {
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
        auto model = molecule.get_histogram();
        auto* exv_model = dynamic_cast<hist::IScatteringModelExv*>(model.get());
        REQUIRE(exv_model != nullptr);
        exv_model->apply_water_scaling_factor(3);
        exv_model->apply_excluded_volume_scaling_factor(1.05);
        exv_model->apply_solvent_density_scaling_factor(1.1);
        exv_model->apply_atomic_debye_waller_factor(0.5);
        exv_model->apply_exv_debye_waller_factor(0.3);
        auto I = model->debye_transform();
        auto q = I.get_axis().as_vector();
        auto A = molecule.get_scattering_amplitude()->evaluate(q, std::nullopt, model->get_parameters());
        check_consistent(q, A, I.get_intensity());

        // and it differs from the amplitude at the default parameters
        auto A_default = molecule.get_scattering_amplitude()->evaluate(q);
        REQUIRE_THAT(A.front(), !Catch::Matchers::WithinRel(A_default.front(), 1e-2));
    }
}

TEST_CASE("Molecule::get_scattering_amplitude: invariant under translation") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
    Molecule molecule("tests/files/2epe.pdb");
    molecule.generate_new_hydration();

    std::vector<double> q = {0.01, 0.1, 0.3, 0.5};
    auto A1 = molecule.get_scattering_amplitude()->evaluate(q);
    molecule.translate({10, -5, 3});
    auto A2 = molecule.get_scattering_amplitude()->evaluate(q);
    for (std::size_t i = 0; i < q.size(); ++i) {
        REQUIRE_THAT(A2[i], Catch::Matchers::WithinRel(A1[i], 1e-4));
    }
}

// moving the centre a distance d from the centroid lowers the amplitude by the factor 1 - q^2 d^2/6 at small q
TEST_CASE("Molecule::get_scattering_amplitude: the centre can be chosen") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    Molecule molecule("tests/files/2epe.pdb");
    molecule.generate_new_hydration();

    // the default centre lies where the q^2 term is smallest, so a shift in any direction lowers the amplitude
    std::vector<double> q = {0.002, 0.005};
    auto A = molecule.get_scattering_amplitude()->evaluate(q);
    auto centroid = molecule.get_scattering_amplitude()->default_centre();
    auto A_centroid = molecule.get_scattering_amplitude()->evaluate(q, centroid);
    for (std::size_t i = 0; i < q.size(); ++i) {REQUIRE_THAT(A_centroid[i], Catch::Matchers::WithinRel(A[i], 1e-9));}

    constexpr double d = 10;
    auto A_shifted = molecule.get_scattering_amplitude()->evaluate(q, centroid + Vector3<double>(0, d, 0));
    for (std::size_t i = 0; i < q.size(); ++i) {
        INFO("q = " << q[i]);
        REQUIRE_THAT(A_shifted[i]/A[i], Catch::Matchers::WithinAbs(1 - q[i]*q[i]*d*d/6, 5e-6));
    }
}

TEST_CASE("Molecule::get_scattering_amplitude: unsupported excluded volume models throw") {
    settings::general::verbose = false;
    settings::exv::exv_method = GENERATE(settings::exv::ExvMethod::Grid, settings::exv::ExvMethod::Average, settings::exv::ExvMethod::CRYSOL);
    Molecule molecule("tests/files/2epe.pdb");
    REQUIRE_THROWS(molecule.get_scattering_amplitude());
}
