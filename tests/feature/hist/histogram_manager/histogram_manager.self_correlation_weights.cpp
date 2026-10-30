#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <settings/All.h>

#include <memory>

using namespace ausaxs;
using namespace ausaxs::data;
using ff_t = form_factor::form_factor_t;

namespace {
    using FFHistogram = hist::CompositeDistanceHistogramFFAvgBase<form_factor::lookup::table_t>;

    // two coincident C atoms and two coincident waters, whose distinct pairs share the first bin with the self-correlations
    Molecule make_molecule() {
        return Molecule({
            Body{
                std::vector{AtomFF({0, 0, 0}, ff_t::C), AtomFF({0, 0, 0}, ff_t::C), AtomFF({3, 0, 0}, ff_t::O)},
                std::vector{Water({0, 3, 0}), Water({0, 3, 0})}
            },
            Body{std::vector{AtomFF({10, 0, 0}, ff_t::C)}, std::vector{Water({10, 3, 0})}}
        });
    }

    struct Expected {
        double self_c, bin0_c;   // the self-correlation weight of the C slot, and the first bin of its diagonal partial
        double self_o, bin0_o;   // the same for the O slot
        double self_ww, bin0_ww; // the same for the waters
    };

    void check(std::unique_ptr<hist::ICompositeDistanceHistogram> result, const Expected& expected) {
        auto* ff = dynamic_cast<FFHistogram*>(result.get());
        REQUIRE(ff != nullptr);
        auto map = form_factor::manager::get_active_mapping();
        int c = map[static_cast<int>(ff_t::C)], o = map[static_cast<int>(ff_t::O)];

        const auto& p_aa = ff->get_raw_aa_counts_by_ff();
        CHECK(ff->get_aa_self_weights()[c] == expected.self_c);
        CHECK(p_aa.index(c, c, 0) == expected.bin0_c);
        CHECK(ff->get_aa_self_weights()[o] == expected.self_o);
        CHECK(p_aa.index(o, o, 0) == expected.bin0_o);
        CHECK(ff->get_ww_self_weight() == expected.self_ww);
        CHECK(ff->get_raw_ww_counts_by_ff().index(0) == expected.bin0_ww);
    }
}

TEST_CASE("HistogramManager: self-correlation weights exclude coincident distinct points") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;

    auto molecule = make_molecule();
    form_factor::manager::use_form_factors(molecule); // the managers are constructed directly, bypassing the factory which normally selects these

    // 3 C atoms, of which 2 coincide; 1 O atom; 3 waters, of which 2 coincide. a coincident pair counts twice, once in either order
    Expected expected{.self_c = 3, .bin0_c = 3+2, .self_o = 1, .bin0_o = 1, .self_ww = 3, .bin0_ww = 3+2};

    SECTION("HistogramManagerMTFFAvg") {
        check(hist::HistogramManagerMTFFAvg<true>(&molecule).calculate_all(), expected);
        check(hist::HistogramManagerMTFFAvg<false>(&molecule).calculate_all(), expected);
    }

    SECTION("HistogramManagerMTFFExplicit") {
        check(hist::HistogramManagerMTFFExplicit<true>(&molecule).calculate_all(), expected);
    }

    SECTION("HistogramManagerMTFFGrid") {
        check(hist::HistogramManagerMTFFGrid(&molecule).calculate_all(), expected);
    }

    SECTION("PartialHistogramManagerMTFF") {
        hist::PartialHistogramManagerMTFF<true> manager(&molecule);
        check(manager.calculate_all(), expected);
        check(manager.calculate_all(), expected); // the self-correlations are not recalculated, so their weights must persist
    }

    SECTION("symmetry managers") {
        // one copy of the second body, far away from everything else
        molecule.get_body(1).symmetry().add(std::make_unique<symmetry::CyclicSymmetry>(
            Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 50}, Vector3<double>{0, 0, 1}, 0, 1
        ));
        expected.self_c += 1;
        expected.bin0_c += 1;

        check(hist::SymmetryManagerMTFF<true>(&molecule).calculate_all(), expected);

        hist::PartialSymmetryManagerMTFF<true> manager(&molecule);
        check(manager.calculate_all(), expected);
        check(manager.calculate_all(), expected);
    }
}
