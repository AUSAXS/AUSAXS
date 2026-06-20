#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "settings/HistogramSettings.h"

#include <random>

using namespace ausaxs;
using namespace ausaxs::data;

auto test_reference_symmetry = [] (settings::hist::HistogramManagerChoice choice) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> d(-8, 8);

    // a cyclic symmetry shared by two bodies, replicating the {body0, body1} pair as a unit
    auto[angle, reps] = GENERATE(
        std::make_pair(std::numbers::pi, 1),         // shared c2
        std::make_pair(2*std::numbers::pi/3, 2)      // shared c3
    );
    int n_atoms = GENERATE(1, 3);

    for (int i = 0; i < 3; ++i) {
        auto random_atoms = [&](int n) {
            std::vector<AtomFF> atoms;
            for (int j = 0; j < n; ++j) {atoms.push_back(AtomFF({d(gen), d(gen), d(gen)}, form_factor::form_factor_t::C));}
            return atoms;
        };
        Molecule m({Body{random_atoms(n_atoms)}, Body{random_atoms(n_atoms+1)}});
        m.set_histogram_manager(choice);
        set_unity_charge(m);

        symmetry::CyclicSymmetry base(
            symmetry::CyclicSymmetry::_Relation{{6, 0, 0}},
            symmetry::CyclicSymmetry::_Repeat{{0, 0, 1}, angle},
            reps
        );
        // body 0 owns the shared symmetry; body 1 holds a view onto it (located by body+slot)
        m.get_body(0).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(base, std::vector<int>{0, 1}, std::vector<int>{0, 0}, &m));
        m.get_body(1).symmetry().add(std::make_unique<symmetry::ReferenceSymmetryView>(&m, 0, 0));

        auto h = m.get_histogram()->get_weighted_counts();

        // ground truth: materialise every copy of both bodies explicitly
        auto b0 = m.get_body(0).symmetry().explicit_structure();
        auto b1 = m.get_body(1).symmetry().explicit_structure();
        Molecule m2({Body{std::move(b0.atoms), std::move(b0.waters)}, Body{std::move(b1.atoms), std::move(b1.waters)}});
        set_unity_charge(m2);
        auto h2 = m2.get_histogram()->get_weighted_counts();

        CHECK(compare_hist_approx(h, h2));
    }
};

TEST_CASE("SymmetryManager: ReferenceSymmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        test_reference_symmetry(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        test_reference_symmetry(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}