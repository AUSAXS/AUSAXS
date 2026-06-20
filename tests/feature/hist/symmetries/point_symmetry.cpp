#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/TetrahedralSymmetry.h>
#include <data/symmetry/OctahedralSymmetry.h>
#include <data/symmetry/IcosahedralSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "settings/HistogramSettings.h"

using namespace ausaxs;
using namespace ausaxs::data;

auto make_unique_point_sym = [] (const Vector3<double>& translation, const Vector3<double>& rotation) {
    return std::make_unique<symmetry::PointSymmetry>(translation, rotation);
};

auto test_point_symmetry = [] (settings::hist::HistogramManagerChoice choice) {
    SECTION("one body with one atom") {
        // Atom at its own cm so copy is always at atom_pos + translation, independent of rotation
        AtomFF a({0, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a}}});
        m.set_histogram_manager(choice);
        set_unity_charge(m);

        SECTION("single copy via translation") {
            m.get_body(0).symmetry().add(make_unique_point_sym({1, 0, 0}, {0, 0, 0}));
            auto h = m.get_histogram()->get_weighted_counts();
            check_hist(h, {RES(0, 2), RES(1, 2)});
        }

        SECTION("two copies via translation") {
            m.get_body(0).symmetry().add(make_unique_point_sym({1, 0, 0}, {0, 0, 0}));
            m.get_body(0).symmetry().add(make_unique_point_sym({2, 0, 0}, {0, 0, 0}));
            // atoms: {0,0,0}, {1,0,0}, {2,0,0}
            auto h = m.get_histogram()->get_weighted_counts();
            check_hist(h, {RES(0, 3), RES(1, 4), RES(2, 2)});
        }

        SECTION("three copies via translation") {
            m.get_body(0).symmetry().add(make_unique_point_sym({1, 0, 0}, {0, 0, 0}));
            m.get_body(0).symmetry().add(make_unique_point_sym({2, 0, 0}, {0, 0, 0}));
            m.get_body(0).symmetry().add(make_unique_point_sym({3, 0, 0}, {0, 0, 0}));
            // atoms: {0,0,0}, {1,0,0}, {2,0,0}, {3,0,0}
            auto h = m.get_histogram()->get_weighted_counts();
            check_hist(h, {RES(0, 4), RES(1, 6), RES(2, 4), RES(3, 2)});
        }

        SECTION("rotation only (atom at cm: copy overlaps original)") {
            // rotation={0,0,1} → Rz(1 rad); atom at cm={0,0,0}, no translation → copy lands on original
            // 2 atoms at same position: 2 self-counts + 2 cross-pair counts = 4 at bin 0
            m.get_body(0).symmetry().add(make_unique_point_sym({0, 0, 0}, {0, 0, 1}));
            auto h = m.get_histogram()->get_weighted_counts();
            check_hist(h, {RES(0, 4)});
        }

        SECTION("combined translation and rotation (explicit_structure comparison)") {
            m.get_body(0).symmetry().add(make_unique_point_sym({2, 1, 0}, {0, 1, 0}));
            auto h = m.get_histogram()->get_weighted_counts();

            auto b = m.get_body(0).symmetry().explicit_structure();
            auto m2 = Molecule({Body{std::move(b.atoms), std::move(b.waters)}});
            set_unity_charge(m2);
            auto h2 = m2.get_histogram()->get_weighted_counts();
            CHECK(compare_hist_approx(h, h2));
        }
    }

    SECTION("two bodies with one atom each") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({1, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a1}}, Body{std::vector{a2}}});
        m.set_histogram_manager(choice);
        set_unity_charge(m);

        SECTION("single copy of body1") {
            m.get_body(0).symmetry().add(make_unique_point_sym({0, 1, 0}, {0, 0, 0}));
            // atoms: {0,0,0}, {1,0,0}, copy: {0,1,0}
            auto h = m.get_histogram()->get_weighted_counts();
            check_hist(h, {RES(0, 3), RES(1, 4), RES(std::sqrt(2.0), 2)});
        }

        SECTION("copies on both bodies") {
            m.get_body(0).symmetry().add(make_unique_point_sym({0, 1, 0}, {0, 0, 0}));
            m.get_body(1).symmetry().add(make_unique_point_sym({0, 1, 0}, {0, 0, 0}));
            // atoms: {0,0,0}, {1,0,0}, copies: {0,1,0}, {1,1,0}
            auto h = m.get_histogram()->get_weighted_counts();
            check_hist(h, {
                RES(0, 4),
                RES(1, 8),
                RES(std::sqrt(2.0), 4)
            });
        }

        SECTION("rotation with explicit_structure comparison") {
            m.get_body(0).symmetry().add(make_unique_point_sym({3, 0, 0}, {0, 0, 1}));
            auto h = m.get_histogram()->get_weighted_counts();

            auto b = m.get_body(0).symmetry().explicit_structure();
            Molecule m2({Body{std::move(b.atoms), std::move(b.waters)}, Body{std::vector{a2}}});
            set_unity_charge(m2);
            auto h2 = m2.get_histogram()->get_weighted_counts();
            CHECK(compare_hist_approx(h, h2));
        }
    }

    SECTION("multi-atom body with rotation (explicit_structure comparison)") {
        // Non-trivial: cm is between atoms so rotation moves copies away from simple translation
        AtomFF a1({ 1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({-1, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a1, a2}}});
        m.set_histogram_manager(choice);
        set_unity_charge(m);

        SECTION("rotation about y-axis") {
            m.get_body(0).symmetry().add(make_unique_point_sym({5, 0, 0}, {0, 1, 0}));
            auto h = m.get_histogram()->get_weighted_counts();

            auto b = m.get_body(0).symmetry().explicit_structure();
            auto m2 = Molecule({Body{std::move(b.atoms), std::move(b.waters)}});
            set_unity_charge(m2);
            auto h2 = m2.get_histogram()->get_weighted_counts();
            CHECK(compare_hist_approx(h, h2));
        }

        SECTION("rotation about z-axis") {
            m.get_body(0).symmetry().add(make_unique_point_sym({0, 3, 0}, {0, 0, 1}));
            auto h = m.get_histogram()->get_weighted_counts();

            auto b = m.get_body(0).symmetry().explicit_structure();
            auto m2 = Molecule({Body{std::move(b.atoms), std::move(b.waters)}});
            set_unity_charge(m2);
            auto h2 = m2.get_histogram()->get_weighted_counts();
            CHECK(compare_hist_approx(h, h2));
        }
    }
};

TEST_CASE("SymmetryManager: PointSymmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        test_point_symmetry(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        test_point_symmetry(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}