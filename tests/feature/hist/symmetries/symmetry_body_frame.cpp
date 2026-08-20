#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::data;

// A symmetry is anchored to the body it belongs to, so any rigid motion of the base body must carry its copies along
// with it. The assembly - and therefore the scattering - is then invariant under real-space moves, leaving the symmetry
// parameters as the only ones that can change it.
namespace {
    std::vector<AtomFF> make_atoms(const Vector3<double>& offset = {0, 0, 0}) {
        std::vector<AtomFF> atoms = {
            AtomFF({ 1.0,  0.0,  0.0}, form_factor::form_factor_t::C),
            AtomFF({-1.0,  0.5,  0.2}, form_factor::form_factor_t::C),
            AtomFF({ 0.3, -1.4,  0.9}, form_factor::form_factor_t::C),
            AtomFF({ 2.1,  0.7, -1.1}, form_factor::form_factor_t::C)
        };
        for (auto& a : atoms) {a.coordinates() += offset;}
        return atoms;
    }

    // an arbitrary rigid motion: a rotation about an arbitrary pivot, followed by a translation. Whoever moves a body owns its orientation,
    // so the new one is reported to the symmetries exactly as TransformStrategy::restore_symmetry does during a refinement.
    void move_rigidly(Body& b, const Vector3<double>& pivot) {
        auto R = matrix::rotation_matrix<double>({0.4, -0.9, 0.3});
        b.translate(-pivot);
        b.rotate(R);
        b.translate(pivot + Vector3<double>{3, -2, 5});
        b.symmetry().set_orientation(R);
    }

    // two bodies sharing a single c3 about their combined centre: body 0 owns the symmetry, body 1 holds a view of it
    Molecule make_reference_molecule();

    Molecule make_molecule(std::vector<Body>&& bodies) {
        Molecule m(std::move(bodies));
        m.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
        set_unity_charge(m);
        return m;
    }

    Molecule make_reference_molecule() {
        auto m = make_molecule({Body{make_atoms()}, Body{make_atoms({5, 1, -2})}});
        auto base = std::make_unique<symmetry::CyclicSymmetry>(
            Vector3<double>{20, 0, 0}, Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 1}, 2*std::numbers::pi/3, 2
        );
        m.get_body(0).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(std::move(base), std::vector{0, 1}, std::vector{0, 0}, &m));
        m.get_body(1).symmetry().add(std::make_unique<symmetry::ReferenceSymmetryView>(&m, 0, 0));
        return m;
    }
}

TEST_CASE("Symmetry: rigid motion of the base body leaves the assembly unchanged") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    SECTION("PointSymmetry") {
        auto m = make_molecule({Body{make_atoms()}});
        m.get_body(0).symmetry().add(std::make_unique<symmetry::PointSymmetry>(Vector3<double>{8, 0, 0}, Vector3<double>{0.3, 0.7, 0.2}));

        auto h1 = m.get_histogram()->get_weighted_counts();
        move_rigidly(m.get_body(0), m.get_body(0).get_cm());
        CHECK(compare_hist_approx(h1, m.get_histogram()->get_weighted_counts()));
    }

    SECTION("PointSymmetry, rotated about a pivot away from the cm") {
        auto m = make_molecule({Body{make_atoms()}});
        m.get_body(0).symmetry().add(std::make_unique<symmetry::PointSymmetry>(Vector3<double>{8, 0, 0}, Vector3<double>{0.3, 0.7, 0.2}));

        auto h1 = m.get_histogram()->get_weighted_counts();
        move_rigidly(m.get_body(0), {-11, 4, 7});
        CHECK(compare_hist_approx(h1, m.get_histogram()->get_weighted_counts()));
    }

    SECTION("CyclicSymmetry") {
        auto m = make_molecule({Body{make_atoms()}});
        m.get_body(0).symmetry().add(std::make_unique<symmetry::CyclicSymmetry>(
            Vector3<double>{8, 0, 0}, Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 1}, 2*std::numbers::pi/3, 2
        ));

        auto h1 = m.get_histogram()->get_weighted_counts();
        move_rigidly(m.get_body(0), m.get_body(0).get_cm());
        CHECK(compare_hist_approx(h1, m.get_histogram()->get_weighted_counts()));
    }

    SECTION("polyhedral and dihedral symmetries") {
        auto name = GENERATE(std::string("tetrahedral"), "octahedral", "icosahedral", "d3", "dp4");
        auto m = make_molecule({Body{make_atoms()}});
        auto sym = symmetry::create(name);
        symmetry::for_each_leaf(*sym, [] (symmetry::ISymmetry& leaf) {
            auto t = leaf.span_translation();
            auto r = leaf.span_rotation();
            t[0] = 9;
            r[0] = 0.6;
        });
        m.get_body(0).symmetry().add(std::move(sym));

        auto h1 = m.get_histogram()->get_weighted_counts();
        move_rigidly(m.get_body(0), m.get_body(0).get_cm());
        CHECK(compare_hist_approx(h1, m.get_histogram()->get_weighted_counts()));
    }

    SECTION("CompositeSymmetry") {
        auto m = make_molecule({Body{make_atoms()}});
        m.get_body(0).symmetry().add(std::make_unique<symmetry::CompositeSymmetry>(
            std::make_unique<symmetry::PointSymmetry>(Vector3<double>{6, 0, 0}, Vector3<double>{0.3, 0.7, 0.2}),
            std::make_unique<symmetry::CyclicSymmetry>(
                Vector3<double>{15, 0, 0}, Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 1}, 2*std::numbers::pi/3, 2
            )
        ));

        auto h1 = m.get_histogram()->get_weighted_counts();
        move_rigidly(m.get_body(0), m.get_body(0).get_cm());
        CHECK(compare_hist_approx(h1, m.get_histogram()->get_weighted_counts()));
    }

    SECTION("ReferenceSymmetry: the whole group moves rigidly") {
        auto m = make_reference_molecule();

        auto h1 = m.get_histogram()->get_weighted_counts();
        Vector3<double> pivot = {-3, 2, 1};
        move_rigidly(m.get_body(0), pivot);
        move_rigidly(m.get_body(1), pivot);
        CHECK(compare_hist_approx(h1, m.get_histogram()->get_weighted_counts()));
    }
}

// A shared symmetry replicates the participating bodies as one asymmetric unit, so every copy must be placed by a single operator acting on all of them.
// If each body reoriented its own copies, the participants of a copy would drift apart as soon as their orientations diverged - which they do, since the
// optimiser moves each body independently - and the assembly would stop being symmetric.
TEST_CASE("ReferenceSymmetry: copies stay congruent when participants are moved independently") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    auto m = make_reference_molecule();

    // only body 0 is reoriented, so the two participants no longer share an orientation
    auto R = matrix::rotation_matrix<double>({0.4, -0.9, 0.3});
    auto pivot = m.get_body(0).get_cm();
    m.get_body(0).translate(-pivot);
    m.get_body(0).rotate(R);
    m.get_body(0).translate(pivot);
    m.get_body(0).symmetry().set_orientation(R);

    auto A = m.get_body(0).symmetry().explicit_structure(); // [body | copy 1 | copy 2], in that order
    auto B = m.get_body(1).symmetry().explicit_structure();
    std::size_t nA = m.get_body(0).size_atom(), nB = m.get_body(1).size_atom();

    // copy k of body 0 and copy k of body 1 must be related to the two bodies by the same isometry, so every cross distance is preserved copy-for-copy
    for (std::size_t k = 1; k <= 2; ++k) {
        for (std::size_t i = 0; i < nA; ++i) {
            for (std::size_t j = 0; j < nB; ++j) {
                double d0 = A.atoms[i].coordinates().distance(B.atoms[j].coordinates());
                double dk = A.atoms[k*nA+i].coordinates().distance(B.atoms[k*nB+j].coordinates());
                INFO("copy " << k << ", atoms " << i << " and " << j);
                CHECK_THAT(dk, Catch::Matchers::WithinAbs(d0, 1e-9));
            }
        }
    }
}
