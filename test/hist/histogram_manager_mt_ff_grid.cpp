#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <settings/All.h>
#include <utility/Utility.h>

#include "hist_test_helper.h"

using namespace hist;
using namespace data;
using namespace data::record;

template<bool weighted>
auto test = [] (const Molecule& protein) {
    settings::molecule::center = false; // to avoid rounding errors
    auto h = hist::HistogramManagerMTFFGrid<weighted>(&protein).calculate_all();

    // convert the grid to water atoms
    auto exv_grid = protein.get_grid()->generate_excluded_volume();
    std::vector<Water> waters(exv_grid.size());
    for (unsigned int i = 0; i < exv_grid.size(); i++) {
        waters[i] = Water(i, "C", "", "LYS", 'A', 1, "", exv_grid[i], 1, 0, constants::atom_t::dummy, "");
        waters[i].set_effective_charge(1);
    }
    Molecule exv(protein.get_atoms(), waters);
    auto h_exv = hist::HistogramManager<weighted>(&exv).calculate_all();

    // calculate the xx, ax, aa distributions
    auto h_cast = static_cast<CompositeDistanceHistogramFFAvg*>(h.get());
    hist::Distribution1D xx1, ax1, aa1;
    {
        // we have to sum all form factor contributions into the single distributions
        auto aa = h_cast->get_aa_counts_ff();
        hist::Distribution1D temp_aa(aa.size_z()), temp_ax(aa.size_z()), temp_xx(aa.size_z());
        for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < form_factor::get_count_without_excluded_volume(); ++j) {
                std::transform(aa.begin(i, j), aa.end(i, j), temp_aa.begin(), temp_aa.begin(), std::plus<>());
            }
            std::transform(aa.begin(i, form_factor::exv_bin), aa.end(i, form_factor::exv_bin), temp_ax.begin(), temp_ax.begin(), std::plus<>());
            std::transform(aa.begin(form_factor::exv_bin, i), aa.end(form_factor::exv_bin, i), temp_ax.begin(), temp_ax.begin(), std::plus<>()); // should all be 0
        }
        std::transform(aa.begin(form_factor::exv_bin, form_factor::exv_bin), aa.end(form_factor::exv_bin, form_factor::exv_bin), temp_xx.begin(), temp_xx.begin(), std::plus<>());

        aa1 = temp_aa;
        ax1 = temp_ax;
        xx1 = temp_xx;
    }

    auto xx2 = h_exv->get_ww_counts();
    auto ax2 = h_exv->get_aw_counts();
    auto aa2 = h_exv->get_aa_counts();

    if (xx1.size() < xx2.size()) {
        for (unsigned int i = xx1.size(); i < xx2.size(); ++i) {
            REQUIRE(xx2.index(i) == 0);
            REQUIRE(ax2.index(i) == 0);
            REQUIRE(aa2.index(i) == 0);
        }
    } else {
        for (unsigned int i = xx2.size(); i < xx1.size(); ++i) {
            REQUIRE(xx1.index(i) == 0);
            REQUIRE(ax1.index(i) == 0);
            REQUIRE(aa1.index(i) == 0);
        }
    }

    // aa
    for (int i = 0; i < std::min<int>(xx1.size(), xx2.size()); ++i) {
        if (!utility::approx(aa1.index(i), aa2.index(i), 1e-6, 0)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE_THAT(aa1.index(i), Catch::Matchers::WithinAbs(aa2.index(i), 1e-6));
        }
    }
    SUCCEED();

    // xx
    for (int i = 0; i < std::min<int>(xx1.size(), xx2.size()); ++i) {
        if (!utility::approx(xx1.index(i), xx2.index(i), 1e-6, 0)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE_THAT(xx1.index(i), Catch::Matchers::WithinAbs(xx2.index(i), 1e-6));
        }
    } 
    SUCCEED();

    // ax
    for (int i = 0; i < std::min<int>(xx1.size(), xx2.size()); ++i) {
        if (!utility::approx(ax1.index(i), ax2.index(i), 1e-6, 0)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE_THAT(ax1.index(i), Catch::Matchers::WithinAbs(ax2.index(i), 1e-6));
        }
    }
    SUCCEED();
};

// Check that the histograms are correct
TEST_CASE("HistogramManagerMTFFGrid::calculate") {
    settings::molecule::use_effective_charge = false;
    SECTION("simple") {
        settings::grid::width = GENERATE(0.2, 0.5, 1, 2);
        settings::grid::exv_radius = settings::grid::width;

        SECTION(std::string("width = ") + std::to_string(settings::grid::width)) {
            Atom a1(0, "C", "", "LYS", 'A', 1, "", {0, 0, 0}, 1, 0, constants::atom_t::dummy, "");
            Molecule protein({a1});
            test<false>(protein);
            test<true>(protein);
        }
    }

    SECTION("actual data") {
        settings::general::verbose = false;
        settings::grid::width = 1;
        settings::grid::exv_radius = 1;
        Molecule protein("test/files/LAR1-2.pdb");
        protein.clear_hydration();
        test<false>(protein);
        test<true>(protein);
    }
}

// Check that the weighted bins are correct and separate for the excluded volume and the protein atoms
TEST_CASE("HistogramManagerMTFFGrid: weighted_bins", "[files]") {
    settings::molecule::use_effective_charge = false;
    settings::hist::weighted_bins = true;
    Molecule protein("test/files/2epe.pdb");

    SECTION("simple") {
        std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Molecule protein(a);
        set_unity_charge(protein);

        auto exv_grid = protein.get_grid()->generate_excluded_volume();
        std::vector<Atom> atoms(exv_grid.size());
        for (unsigned int i = 0; i < exv_grid.size(); i++) {
            atoms[i] = Atom(i, "C", "", "LYS", 'A', 1, "", exv_grid[i], 1, 0, constants::atom_t::dummy, "");
        }
        Molecule exv(atoms);

        auto h_grid = hist::HistogramManagerMTFFGrid<true>(&protein).calculate_all();
        auto h_exv =  hist::HistogramManagerMT<true>(&exv).calculate_all();
        auto h_atom = hist::HistogramManagerMT<true>(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        CHECK(compare_hist(h_grid_cast->get_d_axis(),   h_atom->get_d_axis()));
        CHECK(compare_hist(h_grid_cast->get_d_axis_xx(), h_exv->get_d_axis()));
    }

    SECTION("real data") {
        auto exv_grid = protein.get_grid()->generate_excluded_volume();
        std::vector<Atom> atoms(exv_grid.size());
        for (unsigned int i = 0; i < exv_grid.size(); i++) {
            atoms[i] = Atom(i, "C", "", "LYS", 'A', 1, "", exv_grid[i], 1, 0, constants::atom_t::dummy, "");
        }
        Molecule exv(atoms);
        for (auto& b : exv.get_bodies()) {for (auto& a : b.get_atoms()) {a.set_effective_charge(1);}}

        auto h_grid = hist::HistogramManagerMTFFGrid<true>(&protein).calculate_all();
        auto h_exv =  hist::HistogramManagerMT<true>(&exv).calculate_all();
        auto h_atom = hist::HistogramManagerMT<true>(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        CHECK(compare_hist(h_grid_cast->get_d_axis(),   h_atom->get_d_axis()));
        CHECK(compare_hist(h_grid_cast->get_d_axis_xx(), h_exv->get_d_axis()));
    }
}