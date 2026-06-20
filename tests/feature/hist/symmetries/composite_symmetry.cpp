#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
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

auto test_composite_symmetry = [] (settings::hist::HistogramManagerChoice choice) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> d(-8, 8);

    // build inner/outer sub-symmetries; the nesting "inner-outer" replicates the inner unit
    enum class kind {p2_c3, c2_c3, c3_c2};
    auto nesting = GENERATE(kind::p2_c3, kind::c2_c3, kind::c3_c2);
    int n_atoms = GENERATE(1, 3);

    auto make_cyclic = [](double angle, int reps, Vector3<double> offset) {
        return std::make_unique<symmetry::CyclicSymmetry>(
            symmetry::CyclicSymmetry::_Relation{offset},
            symmetry::CyclicSymmetry::_Repeat{{0, 0, 1}, angle},
            reps
        );
    };

    for (int i = 0; i < 3; ++i) {
        std::vector<AtomFF> atoms;
        for (int j = 0; j < n_atoms; ++j) {
            atoms.push_back(AtomFF({d(gen), d(gen), d(gen)}, form_factor::form_factor_t::C));
        }

        Molecule m({Body{atoms}});
        m.set_histogram_manager(choice);
        set_unity_charge(m);

        std::unique_ptr<symmetry::ISymmetry> inner, outer;
        switch (nesting) {
            case kind::p2_c3:
                inner = std::make_unique<symmetry::PointSymmetry>(Vector3<double>{4, 1, 0}, Vector3<double>{0, 0, 0});
                outer = make_cyclic(2*std::numbers::pi/3, 2, {6, 0, 0});
                break;
            case kind::c2_c3:
                inner = make_cyclic(std::numbers::pi, 1, {3, 0, 0});
                outer = make_cyclic(2*std::numbers::pi/3, 2, {7, 0, 0});
                break;
            case kind::c3_c2:
                inner = make_cyclic(2*std::numbers::pi/3, 2, {3, 0, 0});
                outer = make_cyclic(std::numbers::pi, 1, {8, 0, 0});
                break;
        }
        m.get_body(0).symmetry().add(std::make_unique<symmetry::CompositeSymmetry>(std::move(inner), std::move(outer)));

        auto h = m.get_histogram()->get_weighted_counts();

        auto b = m.get_body(0).symmetry().explicit_structure();
        auto m2 = Molecule({Body{std::move(b.atoms), std::move(b.waters)}});
        set_unity_charge(m2);
        auto h2 = m2.get_histogram()->get_weighted_counts();

        CHECK(compare_hist_approx(h, h2));
    }
};

TEST_CASE("SymmetryManager: CompositeSymmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        test_composite_symmetry(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        test_composite_symmetry(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}