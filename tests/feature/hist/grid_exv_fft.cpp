#if defined(POCKETFFT_AVAILABLE)

#include <algorithm>
#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <form_factor/FormFactorType.h>
#include <grid/Grid.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/GridExvFFT.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/All.h>

#include <hist/hist_test_helper.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::hist;

namespace {
    // the distance histogram of unit-weight atoms placed on the given points, evaluated by the regular pair-loop manager
    std::unique_ptr<ICompositeDistanceHistogram> pair_histogram(const std::vector<Vector3<double>>& points) {
        std::vector<AtomFF> atoms(points.size());
        std::ranges::transform(points, atoms.begin(), [] (const Vector3<double>& p) {
            return AtomFF{p, form_factor::form_factor_t::C};
        });
        Molecule molecule({Body(atoms)});
        set_unity_charge(molecule);
        molecule.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
        return molecule.get_histogram();
    }

    // the pair counts of a histogram, leaving out the self-pairs the transform does not count
    std::vector<double> pair_counts(const ICompositeDistanceHistogram& h, int self_pairs) {
        auto counts = h.get_aa_counts().as_vector();
        counts[0] -= self_pairs;
        return counts;
    }

    // a pair crossing a bin edge moves the weighted centres of both bins, but keeps each within its bin
    void compare_axes(const ICompositeDistanceHistogram& expected, const WeightedDistribution1D& actual) {
        REQUIRE(compare_hist(expected.get_d_axis(), actual.get_weighted_axis(), settings::axes::bin_width, 0));
    }

    double total_count(const WeightedDistribution1D& p) {
        double sum = 0;
        for (int i = 0; i < p.size(); ++i) {sum += static_cast<double>(p.index(i).count);}
        return sum;
    }
}

// The transform must reproduce the histogram of the regular pair loop. The pair loop forms its distances in float while
// the transform uses double, so pairs close to a bin edge may land in a neighbouring bin.
TEST_CASE("lattice::self_correlation: matches the pair loop", "[files]") {
    settings::general::verbose = false;
    settings::grid::cell_width = 1;
    settings::grid::exv::width = 1;

    Molecule protein("tests/files/2epe.pdb");
    protein.clear_hydration();
    auto exv = grid::exv::RawGridExv::create(protein.get_grid());
    REQUIRE(exv.spacing == 1);
    REQUIRE(!exv.interior.empty());
    REQUIRE(exv.interior_sites.size() == exv.interior.size());

    int bin_count = hist::detail::required_bin_count<false>(exv.interior);
    double inv_bin_width = hist::detail::WidthController<false>::get_inv_width();
    auto lattice = hist::detail::lattice::self_correlation(exv, inv_bin_width, bin_count);
    int n = static_cast<int>(exv.interior.size());
    CHECK(total_count(lattice) == static_cast<double>(n)*(n-1));

    auto expected = pair_histogram(exv.interior);
    REQUIRE(compare_hist_approx(pair_counts(*expected, n), lattice.get_content()));
    compare_axes(*expected, lattice);
}

TEST_CASE("lattice::correlations: matches the pair loops", "[files]") {
    settings::general::verbose = false;
    settings::grid::cell_width = 1;
    settings::grid::exv::width = 1;
    settings::grid::exv::surface_thickness = 1;

    Molecule protein("tests/files/2epe.pdb");
    protein.clear_hydration();
    auto exv = grid::exv::RawGridWithSurfaceExv::create(protein.get_grid());
    REQUIRE(exv.spacing == 1);
    REQUIRE(!exv.interior.empty());
    REQUIRE(!exv.surface.empty());
    REQUIRE(exv.interior_sites.size() == exv.interior.size());
    REQUIRE(exv.surface_sites.size() == exv.surface.size());

    int bin_count = hist::detail::required_bin_count<false>(exv.interior, exv.surface);
    double inv_bin_width = hist::detail::WidthController<false>::get_inv_width();
    auto lattice = hist::detail::lattice::correlations(exv, inv_bin_width, bin_count);

    int n_i = static_cast<int>(exv.interior.size());
    int n_s = static_cast<int>(exv.surface.size());
    auto interior = pair_histogram(exv.interior);
    auto surface = pair_histogram(exv.surface);
    SECTION("interior") {
        REQUIRE(compare_hist_approx(pair_counts(*interior, n_i), lattice.first.get_content()));
        compare_axes(*interior, lattice.first);
    }
    SECTION("surface") {
        REQUIRE(compare_hist_approx(pair_counts(*surface, n_s), lattice.second.get_content()));
        compare_axes(*surface, lattice.second);
    }

    // the cross term is what the combined set holds beyond the pairs within each set
    SECTION("cross") {
        auto all_points = exv.interior;
        all_points.insert(all_points.end(), exv.surface.begin(), exv.surface.end());
        auto cross = pair_counts(*pair_histogram(all_points), 0);
        auto counts_i = pair_counts(*interior, 0);
        auto counts_s = pair_counts(*surface, 0);
        for (int i = 0; i < static_cast<int>(cross.size()); ++i) {
            cross[i] -= bin_or_zero(counts_i, i) + bin_or_zero(counts_s, i);
        }
        REQUIRE(compare_hist_approx(cross, lattice.cross.get_content()));
    }
}

// A filled 4x4x4 cube, small enough to count by hand: every ordered pair is counted exactly once, and the nearest
// neighbours land in the bin of the lattice spacing.
TEST_CASE("lattice::self_correlation: counts a small cube") {
    settings::general::verbose = false;

    grid::exv::GridExcludedVolume exv;
    exv.spacing = 2;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                exv.interior.emplace_back(2*i, 2*j, 2*k);
                exv.interior_sites.emplace_back(i, j, k);
            }
        }
    }

    auto lattice = hist::detail::lattice::self_correlation(exv, 1, 100);
    CHECK(total_count(lattice) == 64*63);

    // 3*4*4 adjacent pairs along each of the three axes, counted in both directions
    CHECK(lattice.index(2).count == 2*3*(3*4*4));
    CHECK(lattice.index(0).count == 0);
    CHECK(lattice.index(1).count == 0);
}

#endif
