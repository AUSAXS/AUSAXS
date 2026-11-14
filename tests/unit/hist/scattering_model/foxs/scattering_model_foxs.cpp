#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/atoms/AtomFF.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/foxs/ScatteringModelFoXS.h>
#include <settings/ExvSettings.h>
#include <settings/ScatteringSettings.h>

using namespace ausaxs;
using namespace ausaxs::data;
using ff_t = form_factor::form_factor_t;

TEST_CASE("ScatteringModelFoXS: refuses neutrons") {
    settings::exv::exv_method = settings::exv::ExvMethod::FoXS;
    Molecule molecule({Body(std::vector<AtomFF>{
        AtomFF({0, 0, 0}, ff_t::CH3),
        AtomFF({1, 0, 0}, ff_t::O),
    })});

    auto construct = [&] () {
        form_factor::manager::use_form_factors(molecule); // the histogram is constructed directly, bypassing the factory which normally selects these
        int n = form_factor::get_active_count();
        return hist::ScatteringModelFoXS(hist::Distribution3D<hist::Shape::Triangular>(n, n, 1), hist::Distribution2D(n, 1), hist::Distribution1D(1), std::vector<double>(n, 0), 0, hist::Distribution1D(1));
    };

    SECTION("X-ray") {
        settings::scattering::radiation = settings::scattering::Radiation::XRay;
        CHECK_NOTHROW(construct());
    }

    SECTION("neutron") {
        settings::scattering::radiation = settings::scattering::Radiation::Neutron;
        CHECK_THROWS(construct());
    }
}
