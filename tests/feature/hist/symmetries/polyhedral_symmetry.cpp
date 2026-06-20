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

#include <random>

using namespace ausaxs;
using namespace ausaxs::data;

auto test_polyhedral_symmetry = [] (settings::hist::HistogramManagerChoice choice) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> d(-10, 10);
    static std::uniform_real_distribution<> r(-std::numbers::pi, std::numbers::pi);

    auto make_group = GENERATE(
        +[] () -> std::unique_ptr<symmetry::IPolyhedralSymmetry> {return std::make_unique<symmetry::TetrahedralSymmetry>();},
        +[] () -> std::unique_ptr<symmetry::IPolyhedralSymmetry> {return std::make_unique<symmetry::OctahedralSymmetry>();},
        +[] () -> std::unique_ptr<symmetry::IPolyhedralSymmetry> {return std::make_unique<symmetry::IcosahedralSymmetry>();}
    );
    int n_atoms = GENERATE(1, 4);

    // build a body, apply a polyhedral symmetry, and compare the reused-distance histogram
    // against the ground truth obtained by explicitly materialising every copy
    for (int i = 0; i < 3; ++i) {
        std::vector<AtomFF> atoms;
        for (int j = 0; j < n_atoms; ++j) {
            atoms.push_back(AtomFF({d(gen), d(gen), d(gen)}, form_factor::form_factor_t::C));
        }

        Molecule m({Body{atoms}});
        m.set_histogram_manager(choice);
        set_unity_charge(m);

        auto sym = make_group();
        sym->translation = {d(gen), d(gen), d(gen)};
        sym->rotation = {r(gen), r(gen), r(gen)};
        m.get_body(0).symmetry().add(std::move(sym));

        auto h = m.get_histogram()->get_weighted_counts();

        auto b = m.get_body(0).symmetry().explicit_structure();
        auto m2 = Molecule({Body{std::move(b.atoms), std::move(b.waters)}});
        set_unity_charge(m2);
        auto h2 = m2.get_histogram()->get_weighted_counts();

        CHECK(compare_hist_approx(h, h2));
    }
};

TEST_CASE("SymmetryManager: PolyhedralSymmetry") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        test_polyhedral_symmetry(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        test_polyhedral_symmetry(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}

auto test_polyhedral_symmetry_lysozyme = [] (settings::hist::HistogramManagerChoice choice) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> d(-10, 10);
    static std::uniform_real_distribution<> r(-std::numbers::pi, std::numbers::pi);

    Molecule m({Body("tests/files/2epe.pdb")});
    m.set_histogram_manager(choice);
    m.generate_new_hydration();
    set_unity_charge(m);

    auto sym = std::make_unique<symmetry::TetrahedralSymmetry>();
    sym->translation = {d(gen), d(gen), d(gen)};
    sym->rotation = {r(gen), r(gen), r(gen)};
    m.get_body(0).symmetry().add(std::move(sym));
    auto h = m.get_histogram()->get_weighted_counts();

    auto b = m.get_body(0).symmetry().explicit_structure();
    auto m2 = Molecule({Body{std::move(b.atoms), m.get_waters()}});
    set_unity_charge(m2);
    auto h2 = m2.get_histogram()->get_weighted_counts();

    CHECK(compare_hist_approx(h, h2));
};

TEST_CASE("SymmetryManager: PolyhedralSymmetry on hydrated lysozyme") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        test_polyhedral_symmetry_lysozyme(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        test_polyhedral_symmetry_lysozyme(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}