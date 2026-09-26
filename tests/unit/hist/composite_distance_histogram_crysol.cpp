#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/atoms/AtomFF.h>
#include <form_factor/ExvTable.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>
#include <settings/ExvSettings.h>

using namespace ausaxs;
using namespace ausaxs::data;
using ff_t = form_factor::form_factor_t;

namespace {
    Molecule make_molecule() {
        return Molecule({Body(std::vector<AtomFF>{
            AtomFF({0, 0, 0}, ff_t::CH3),
            AtomFF({1, 0, 0}, ff_t::O),
            AtomFF({0, 1, 0}, ff_t::NH2),
        })});
    }

    template<typename T>
    T construct(const Molecule& molecule) {
        constexpr int n = form_factor::total_ff_count;
        return T(hist::Distribution3D<hist::Shape::Triangular>(n, n, 1), hist::Distribution2D(n, 1), hist::Distribution1D(1), hist::Distribution1D(1), &molecule);
    }
}

TEST_CASE("CompositeDistanceHistogramCrysol: always uses the Traube volumes") {
    auto original_set = settings::exv::exv_set.value;
    settings::exv::exv_set = settings::exv::ExvSet::vdw;
    auto molecule = make_molecule();
    double vdw_V = form_factor::ExvTableManager::get_average_displaced_volume(&molecule);

    auto check = [&] (const hist::CompositeDistanceHistogramCrysol& h) {
        CHECK(settings::exv::exv_set == settings::exv::ExvSet::Traube);
        CHECK(*form_factor::ExvTableManager::get_current_exv_table() == constants::exv::Traube);

        // the average volume must be evaluated after switching to the Traube volumes
        double traube_V = (constants::exv::Traube.get(ff_t::CH3) + constants::exv::Traube.get(ff_t::O) + constants::exv::Traube.get(ff_t::NH2))/3;
        REQUIRE(traube_V != vdw_V);
        CHECK_THAT(h.average_displaced_V, Catch::Matchers::WithinRel(traube_V, 1e-12));
    };

    SECTION("CRYSOL") {check(construct<hist::CompositeDistanceHistogramCrysol>(molecule));}
    SECTION("Pepsi")  {check(construct<hist::CompositeDistanceHistogramPepsi>(molecule));}

    settings::exv::exv_set = original_set;
}
