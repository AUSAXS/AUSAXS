#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/DihedralSymmetry.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "settings/HistogramSettings.h"

#include <numbers>
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
        m.get_body(0).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(
            std::make_unique<symmetry::CyclicSymmetry>(base), std::vector<int>{0, 1}, std::vector<int>{0, 0}, &m
        ));
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

auto test_reference_symmetry_dihedral = [] (settings::hist::HistogramManagerChoice choice) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> d(-8, 8);
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

        auto base = std::make_unique<symmetry::DihedralSymmetry<2>>();
        base->translation = {6, 0, 0};
        // body 0 owns the shared symmetry; body 1 holds a view onto it (located by body+slot)
        m.get_body(0).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(
            std::move(base), std::vector<int>{0, 1}, std::vector<int>{0, 0}, &m
        ));
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

TEST_CASE("SymmetryManager: ReferenceSymmetry with dihedral base") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        test_reference_symmetry_dihedral(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        test_reference_symmetry_dihedral(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}

auto test_reference_symmetry_after_transform = [] (settings::hist::HistogramManagerChoice choice) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> d(-8, 8);
    static std::uniform_real_distribution<> d_angle(-std::numbers::pi, std::numbers::pi);

    auto[angle, reps] = GENERATE(
        std::make_pair(std::numbers::pi, 1),         // shared c2
        std::make_pair(2*std::numbers::pi/3, 2)      // shared c3
    );

    auto random_atoms = [&](int n) {
        std::vector<AtomFF> atoms;
        for (int j = 0; j < n; ++j) {atoms.push_back(AtomFF({d(gen), d(gen), d(gen)}, form_factor::form_factor_t::C));}
        return atoms;
    };
    Molecule m({Body{random_atoms(3)}, Body{random_atoms(4)}});
    m.set_histogram_manager(choice);
    set_unity_charge(m);

    symmetry::CyclicSymmetry base(
        symmetry::CyclicSymmetry::_Relation{{6, 0, 0}},
        symmetry::CyclicSymmetry::_Repeat{{0, 0, 1}, angle},
        reps
    );
    m.get_body(0).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(
        std::make_unique<symmetry::CyclicSymmetry>(base), std::vector<int>{0, 1}, std::vector<int>{0, 0}, &m
    ));
    m.get_body(1).symmetry().add(std::make_unique<symmetry::ReferenceSymmetryView>(&m, 0, 0));

    auto ground_truth = [&] {
        auto b0 = m.get_body(0).symmetry().explicit_structure();
        auto b1 = m.get_body(1).symmetry().explicit_structure();
        Molecule m2({Body{std::move(b0.atoms), std::move(b0.waters)}, Body{std::move(b1.atoms), std::move(b1.waters)}});
        set_unity_charge(m2);
        return m2.get_histogram()->get_weighted_counts();
    };

    CHECK(compare_hist_approx(m.get_histogram()->get_weighted_counts(), ground_truth()));

    // rigidly rotate+translate BOTH participating bodies by the same transform, about their shared combined cm
    auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(m.get_body(0).symmetry().get(0));
    Vector3<double> pivot = ref->combined_cm();
    Vector3<double> axis{d(gen), d(gen), d(gen)};
    axis = axis/axis.magnitude();
    auto R = matrix::rotation_matrix<double>(axis, d_angle(gen));
    Vector3<double> t{d(gen), d(gen), d(gen)};
    for (int i = 0; i < 2; ++i) {
        auto& body = m.get_body(i);
        body.translate(-pivot);
        body.rotate(R);
        body.translate(pivot + t);
    }

    CHECK(compare_hist_approx(m.get_histogram()->get_weighted_counts(), ground_truth()));
};

TEST_CASE("SymmetryManager: ReferenceSymmetry stays consistent with ground truth after a further rigid transform") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        test_reference_symmetry_after_transform(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        settings::flags::max_bin_count = constants::axes::d_axis.bins;
        test_reference_symmetry_after_transform(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}

TEST_CASE("ReferenceSymmetry: combined centre of mass is mass-weighted, matching a single body's own centre of mass") {
    settings::flags::max_bin_count = constants::axes::d_axis.bins;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    std::vector<AtomFF> atoms_a{AtomFF({0, 0, 0}, form_factor::form_factor_t::C), AtomFF({1, 0, 0}, form_factor::form_factor_t::C), AtomFF({0, 1, 0}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> atoms_b{
        AtomFF({5, 0, 0}, form_factor::form_factor_t::S), AtomFF({5, 1, 0}, form_factor::form_factor_t::S),
        AtomFF({5, 0, 1}, form_factor::form_factor_t::S), AtomFF({5, 1, 1}, form_factor::form_factor_t::S)
    };

    // "plain": one body holding all the atoms, with a plain (non-shared) symmetry
    std::vector<AtomFF> all_atoms = atoms_a;
    all_atoms.insert(all_atoms.end(), atoms_b.begin(), atoms_b.end());
    Molecule m_plain({Body{all_atoms}});
    m_plain.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    set_unity_charge(m_plain);
    m_plain.get_body(0).symmetry().add(std::make_unique<symmetry::CyclicSymmetry>(
        symmetry::CyclicSymmetry::_Relation{{6, 0, 0}}, symmetry::CyclicSymmetry::_Repeat{{0, 0, 1}, std::numbers::pi}, 1
    ));

    // "reference": the same atoms split across two bodies sharing a ReferenceSymmetry with identical parameters
    Molecule m_ref({Body{atoms_a}, Body{atoms_b}});
    m_ref.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    set_unity_charge(m_ref);
    m_ref.get_body(0).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(
        std::make_unique<symmetry::CyclicSymmetry>(
            symmetry::CyclicSymmetry::_Relation{{6, 0, 0}}, symmetry::CyclicSymmetry::_Repeat{{0, 0, 1}, std::numbers::pi}, 1
        ),
        std::vector<int>{0, 1}, std::vector<int>{0, 0}, &m_ref
    ));
    m_ref.get_body(1).symmetry().add(std::make_unique<symmetry::ReferenceSymmetryView>(&m_ref, 0, 0));

    CHECK(compare_hist_approx(m_plain.get_histogram()->get_weighted_counts(), m_ref.get_histogram()->get_weighted_counts()));

    // translating the whole assembly rigidly must preserve the match
    Vector3<double> t{2, -3, 1};
    m_plain.get_body(0).translate(t);
    m_ref.get_body(0).translate(t);
    m_ref.get_body(1).translate(t);

    CHECK(compare_hist_approx(m_plain.get_histogram()->get_weighted_counts(), m_ref.get_histogram()->get_weighted_counts()));
}