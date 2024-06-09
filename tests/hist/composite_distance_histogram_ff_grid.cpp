#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <dataset/SimpleDataset.h>
#include <plots/All.h>
#include <settings/All.h>

#include "../test/hist/hist_test_helper.h"

using namespace data;
using namespace data::record;

TEST_CASE("CompositeDistanceHistogramFFGrid::volumes", "[manual]") {
    settings::general::verbose = false;
    data::Molecule protein("test/files/2epe.pdb");
    auto V = protein.get_volume_grid();

    std::vector<double> volumes;
    std::vector<double> rxs;
    for (double rx = 0.1; rx < 1; rx += 0.1) {
        settings::grid::width = rx;
        settings::grid::exv_radius = rx;
        protein.clear_grid();
        volumes.push_back(protein.get_volume_grid());
        rxs.push_back(rx);

        // hist::CompositeDistanceHistogramFFGrid::regenerate_table();
        // auto h = hist::HistogramManagerMTFFGrid<false>(&protein).calculate_all();
        // auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGrid*>(h.get());
        // auto ffx = h_cast->get_ff_table().index(form_factor::exv_bin, form_factor::exv_bin);

        // auto gridx = protein.get_grid()->generate_excluded_volume();
        // volumes.push_back(gridx.size()*ffx.evaluate(0)/constants::charge::density::water);
        // rxs.push_back(rx);

        // CHECK_THAT(V, Catch::Matchers::WithinRel(protein.get_volume_grid(), 1e-1));
    }

    SimpleDataset dataset(rxs, volumes);
    plots::PlotDataset::quick_plot(dataset, plots::PlotOptions({{"xlabel", "Grid width [Å]"}, {"ylabel", "Volume [Å³]"}, {"color", style::color::blue}}), "temp/test/hist/composite_distance_histogram_ff_grid/volumes.png");
}

TEST_CASE("CompositeDistanceHistogramFFGridSurface: compare_profiles") {}