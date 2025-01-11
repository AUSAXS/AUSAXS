#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <settings/All.h>
#include <utility/Utility.h>

#include "hist_test_helper.h"
#include "grid/grid_debug.h"

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::data;

auto test = [] (const Molecule& protein, std::function<std::unique_ptr<ICompositeDistanceHistogram>(const Molecule&)> calculate) {
    settings::molecule::center = false; // to avoid rounding errors
    auto h = calculate(protein);

    // convert the grid to water atoms
    auto exv_grid = protein.get_grid()->generate_excluded_volume(false);
    std::vector<Water> waters(exv_grid.interior.size());
    for (unsigned int i = 0; i < exv_grid.interior.size(); i++) {
        waters[i] = Water(exv_grid.interior[i]);
        waters[i].weight() = 1;
    }

    Molecule exv({Body{protein.get_atoms(), waters}});
    auto h_exv = hist::HistogramManager<true>(&exv).calculate_all();

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

        aa1 = std::move(temp_aa);
        ax1 = std::move(temp_ax);
        xx1 = std::move(temp_xx);
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

// Check that the Grid histograms are correct
TEST_CASE("HistogramManagerMTFFGrid::calculate", "[files]") {
    SECTION("simple") {
        settings::grid::cell_width = GENERATE(0.2, 0.5, 1, 2);
        settings::grid::exv::width = settings::grid::cell_width;

        SECTION(std::string("width = ") + std::to_string(settings::grid::cell_width)) {
            AtomFF a1({0, 0, 0}, form_factor::form_factor_t::UNKNOWN);
            Molecule protein({Body{std::vector{a1}}});
            test(protein, [](const Molecule& protein) {return hist::HistogramManagerMTFFGrid(&protein).calculate_all();});
        }
    }

    SECTION("actual data") {
        settings::general::verbose = false;
        settings::grid::cell_width = 1;
        settings::grid::exv::width = 1;
        Molecule protein("tests/files/LAR1-2.pdb");
        protein.clear_hydration();
        test(protein, [](const Molecule& protein) {return hist::HistogramManagerMTFFGrid(&protein).calculate_all();});
    }
}

template<typename H, typename C>
auto test_derived = [] () {
    settings::molecule::center = false;
    settings::general::verbose = false;

    SECTION("simple") {
        settings::grid::cell_width = GENERATE(0.2, 0.5, 1, 2);
        settings::grid::exv::width = settings::grid::cell_width;
        settings::grid::exv::surface_thickness = settings::grid::cell_width;

        SECTION(std::string("width = ") + std::to_string(settings::grid::cell_width)) {
            AtomFF a1({0, 0, 0}, form_factor::form_factor_t::UNKNOWN);
            Molecule protein({Body{std::vector{a1}}});
            test(protein, [](const Molecule& protein) {return H(&protein).calculate_all();});
        }
    }

    SECTION("actual data") {
        settings::general::verbose = false;
        settings::grid::cell_width = 1;
        settings::grid::exv::width = 1;
        settings::grid::exv::surface_thickness = 1;
        Molecule protein("tests/files/LAR1-2.pdb");
        protein.clear_hydration();
        test(protein, [](const Molecule& protein) {return H(&protein).calculate_all();});
    }

    SECTION("compare with Grid") {
        settings::grid::min_bins = 20;

        Molecule protein("tests/files/LAR1-2.pdb");
        protein.generate_new_hydration();

        auto h_grid  = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_grids = H(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        auto aa1 = h_grid_cast->get_aa_counts_ff();
        auto ax1 = h_grid_cast->get_aw_counts_ff();
        auto xx1 = h_grid_cast->get_ww_counts_ff();

        auto h_grids_cast = static_cast<C*>(h_grids.get());
        auto aa2 = h_grids_cast->get_aa_counts_ff();
        auto ax2 = h_grids_cast->get_aw_counts_ff();
        auto xx2 = h_grids_cast->get_ww_counts_ff();

        CHECK(xx1.size() == xx2.size());
        for (unsigned int k = 0; k < xx1.size(); ++k) {
            if (!utility::approx(xx1.index(k), xx2.index(k), 1e-3, 0)) {
                std::cout << "histogram_manager_mt_ff_grid failed at index " << k << std::endl;
                REQUIRE_THAT(xx2.index(k), Catch::Matchers::WithinAbs(xx1.index(k), 1e-3));
            }
            SUCCEED();
        }

        CHECK((ax1.size_x() == ax2.size_x() && ax1.size_y() == ax2.size_y()));
        for (unsigned int k = 0; k < ax2.size_y(); ++k) {
            for (unsigned int j = 0; j < ax2.size_x(); ++j) {
                if (!utility::approx(ax2.index(j, k), ax1.index(j, k), 1e-3, 0)) {
                    std::cout << "histogram_manager_mt_ff_grid failed at index " << j << ", " << k << std::endl;
                    REQUIRE_THAT(ax2.index(j, k), Catch::Matchers::WithinAbs(ax1.index(j, k), 1e-3));
                }
                SUCCEED();
            }
        }

        CHECK((aa1.size_x() == aa2.size_x() && aa1.size_y() == aa2.size_y() && aa1.size_z() == aa2.size_z()));
        for (unsigned int k = 0; k < aa2.size_z(); ++k) {
            for (unsigned int j = 0; j < aa2.size_y(); ++j) {
                for (unsigned int i = 0; i < aa2.size_x(); ++i) {
                    if (!utility::approx(aa2.index(i, j, k), aa1.index(i, j, k), 1e-3, 0)) {
                        std::cout << "histogram_manager_mt_ff_grid failed at index " << i << ", " << j << ", " << k << std::endl;
                        REQUIRE_THAT(aa2.index(i, j, k), Catch::Matchers::WithinAbs(aa1.index(i, j, k), 1e-3));
                    }
                    SUCCEED();
                }
            }
        }
    }
};

// Check that the GridSurface histograms are correct
TEST_CASE("HistogramManagerMTFFGridSurface::calculate", "[files]") {
    test_derived<HistogramManagerMTFFGridSurface, CompositeDistanceHistogramFFGridSurface>();
}

TEST_CASE("HistogramManagerMTFFGridScalableExv::calculate", "[files]") {
    test_derived<HistogramManagerMTFFGridScalableExv, CompositeDistanceHistogramFFGridScalableExv>();
}

// Check that the weighted bins are correct and separate for the excluded volume and the protein atoms
TEST_CASE("HistogramManagerMTFFGrid: weighted_bins", "[files]") {
    settings::molecule::center = false;
    settings::hist::weighted_bins = true;
    settings::general::verbose = false;

    std::string file = GENERATE("tests/files/2epe.pdb", "tests/files/LAR1-2.pdb", "tests/files/diamond.pdb", "tests/files/c60.pdb");
    Molecule protein(file);

    SECTION("simple") {
        std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Molecule protein(a);
        set_unity_charge(protein);

        auto exv_grid = protein.get_grid()->generate_excluded_volume(false);
        std::vector<AtomFF> atoms(exv_grid.interior.size());
        for (unsigned int i = 0; i < exv_grid.interior.size(); i++) {
            atoms[i] = AtomFF(exv_grid.interior[i], form_factor::form_factor_t::UNKNOWN);
        }
        Molecule exv({Body{std::vector{atoms}}});

        auto h_grid  = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_grids = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();
        auto h_gridsx= hist::HistogramManagerMTFFGridScalableExv(&protein).calculate_all();
        auto h_exv   = hist::HistogramManagerMT<true>(&exv).calculate_all();
        auto h_atom  = hist::HistogramManagerMT<true>(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        auto h_grids_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(h_grids.get());
        auto h_gridsx_cast = static_cast<CompositeDistanceHistogramFFGridScalableExv*>(h_gridsx.get());

        CHECK(compare_hist(h_grid_cast->get_d_axis(),   h_atom->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis(),  h_atom->get_d_axis()));
        CHECK(compare_hist(h_gridsx_cast->get_d_axis(), h_atom->get_d_axis()));

        CHECK(compare_hist(h_grid_cast->get_d_axis_xx(), h_exv->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis_xx(), h_exv->get_d_axis()));
        CHECK(compare_hist(h_gridsx_cast->get_d_axis_xx(), h_exv->get_d_axis()));
    }

    SECTION("simple, all") {
        settings::grid::min_exv_radius = 0;
        std::vector<AtomFF> atoms = SimpleCube::get_atoms();
        atoms.push_back(AtomFF({0, 0, 0}, form_factor::form_factor_t::C));
        std::for_each(atoms.begin(), atoms.end(), [](AtomFF& a) {a.weight() = 1;});

        Molecule protein({Body{atoms}});
        GridDebug::generate_debug_grid(protein); // overrides exv generation
        auto h = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h.get());

        // check the distance axes
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));
    }

    SECTION("real data") {
        auto exv_grid = protein.get_grid()->generate_excluded_volume(false);
        std::vector<AtomFF> atoms(exv_grid.interior.size());
        for (unsigned int i = 0; i < exv_grid.interior.size(); i++) {
            atoms[i] = AtomFF(exv_grid.interior[i], form_factor::form_factor_t::UNKNOWN);
        }
        Molecule exv({Body{atoms}});
        for (auto& b : exv.get_bodies()) {for (auto& a : b.get_atoms()) {a.weight() = 1;}}

        auto h_grid  = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_grids = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();
        auto h_gridsx= hist::HistogramManagerMTFFGridScalableExv(&protein).calculate_all();
        auto h_exv   = hist::HistogramManagerMT<true>(&exv).calculate_all();
        auto h_atom  = hist::HistogramManagerMT<true>(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        auto h_grids_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(h_grid.get());
        auto h_gridsx_cast = static_cast<CompositeDistanceHistogramFFGridScalableExv*>(h_grid.get());

        CHECK(compare_hist(h_grid_cast->get_d_axis(),   h_atom->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis(),   h_atom->get_d_axis()));
        CHECK(compare_hist(h_gridsx_cast->get_d_axis(),   h_atom->get_d_axis()));

        CHECK(compare_hist(h_grid_cast->get_d_axis_xx(), h_exv->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis_xx(), h_exv->get_d_axis()));
        CHECK(compare_hist(h_gridsx_cast->get_d_axis_xx(), h_exv->get_d_axis()));
    }
}