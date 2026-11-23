#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/scattering_model/AmplitudeCalculator.h>
#include <hist/scattering_model/IScatteringModel.h>
#include <math/ConstexprMath.h>
#include <settings/All.h>

#include <cmath>

using namespace ausaxs;
using namespace ausaxs::data;

// two point scatterers of weights 1 and 3 a distance d apart: the centre is 3d/4 from the first
TEST_CASE("average_amplitude: agrees with the analytical result for two points") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::None;
    constexpr double d = 10;
    Molecule molecule({Body(std::vector<AtomFF>{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C, 1),
        AtomFF({d, 0, 0}, form_factor::form_factor_t::C, 3)
    })});

    std::vector<double> q = {0, 0.01, 0.1, 0.5, 1, 2, 5, 10};
    auto A = hist::average_amplitude(molecule, q);
    for (std::size_t i = 0; i < q.size(); ++i) {
        double expected = std::exp(-0.5*q[i]*q[i])*(math::sinc(q[i]*0.75*d) + 3*math::sinc(q[i]*0.25*d));
        REQUIRE_THAT(A[i], Catch::Matchers::WithinAbs(expected, 1e-9));
    }
}

// the amplitude is the counterpart of the intensity of the same model, so A(0)^2 = I(0), and beta = A^2/I <= 1 at all q since <A>^2 <= <A^2>
TEST_CASE("average_amplitude: consistent with the intensity") {
    settings::general::verbose = false;
    auto exv = GENERATE(
        settings::exv::ExvMethod::Simple, settings::exv::ExvMethod::None,
        settings::exv::ExvMethod::Vacuum, settings::exv::ExvMethod::Fraser
    );
    settings::exv::exv_method = exv;
    SECTION("exv model " + std::to_string(static_cast<int>(exv))) {
        Molecule molecule("tests/files/2epe.pdb");
        molecule.generate_new_hydration();

        auto I = molecule.get_histogram()->debye_transform();
        auto q = I.get_axis().as_vector();
        auto A = hist::average_amplitude(molecule, q);
        REQUIRE(static_cast<int>(A.size()) == I.size());

        CHECK_THAT(A.front()*A.front()/I[0], Catch::Matchers::WithinRel(1, 1e-3));
        for (int i = 0; i < I.size(); ++i) {
            INFO("q = " << q[i] << ", beta = " << A[i]*A[i]/I[i]);
            REQUIRE(A[i]*A[i] <= 1.01*I[i]);
        }
    }
}

TEST_CASE("average_amplitude: invariant under translation") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
    Molecule molecule("tests/files/2epe.pdb");
    molecule.generate_new_hydration();

    std::vector<double> q = {0.01, 0.1, 0.3, 0.5};
    auto A1 = hist::average_amplitude(molecule, q);
    molecule.translate({10, -5, 3});
    auto A2 = hist::average_amplitude(molecule, q);
    for (std::size_t i = 0; i < q.size(); ++i) {
        REQUIRE_THAT(A2[i], Catch::Matchers::WithinRel(A1[i], 1e-4));
    }
}

TEST_CASE("average_amplitude: unsupported excluded volume models throw") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Grid;
    Molecule molecule("tests/files/2epe.pdb");
    REQUIRE_THROWS(hist::average_amplitude(molecule, {0.1}));
}
