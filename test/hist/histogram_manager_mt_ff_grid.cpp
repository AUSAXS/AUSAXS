#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <hydrate/Grid.h>
#include <settings/MoleculeSettings.h>
#include <settings/GridSettings.h>

using namespace hist;
using namespace data;
using namespace data::record;

auto test = [] (const Molecule& protein) {
    auto h = hist::HistogramManagerMTFFGrid<false>(&protein).calculate_all();

    // convert the grid to water atoms
    auto exv_grid = protein.get_grid()->generate_excluded_volume();
    std::vector<Water> waters(exv_grid.size());
    for (unsigned int i = 0; i < exv_grid.size(); i++) {
        waters[i] = Water(i, "C", "", "LYS", 'A', 1, "", exv_grid[i], 1, 0, constants::atom_t::dummy, "");
        waters[i].set_effective_charge(1);
    }
    Molecule exv(protein.get_atoms(), waters);
    auto h_exv = hist::HistogramManagerMT<false>(&exv).calculate_all();

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
            std::transform(aa.begin(form_factor::exv_bin, i), aa.end(form_factor::exv_bin, i), temp_ax.begin(), temp_ax.begin(), std::plus<>());
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

    for (int i = 0; i < std::min<int>(xx1.size(), xx2.size()); ++i) {
        if (xx1.index(i) != xx2.index(i)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE(xx1.index(i) == xx2.index(i));
        } else if (ax1.index(i) != ax2.index(i)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE(ax1.index(i) == ax2.index(i));
        } else if (aa1.index(i) != aa2.index(i)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE(aa1.index(i) == aa2.index(i));
        }
    }
    SUCCEED();    
};

TEST_CASE("HistogramManagerMTFFGrid::calculate") {
    settings::molecule::use_effective_charge = false;

    SECTION("simple") {
        settings::grid::width = GENERATE(0.2, 0.5, 1, 2);
        settings::grid::exv_radius = settings::grid::width;
        SECTION(std::string("width = ") + std::to_string(settings::grid::width)) {
            settings::grid::width = 2;
            settings::grid::exv_radius = 2;
            Atom a1(0, "C", "", "LYS", 'A', 1, "", {0, 0, 0}, 1, 0, constants::atom_t::dummy, "");
            Molecule protein({a1});
            test(protein);
        }
    }

    SECTION("actual data") {
        settings::grid::width = 1;
        settings::grid::exv_radius = 1;
        Molecule protein("test/files/2epe.pdb");
        protein.clear_hydration();
        test(protein);
    }
}