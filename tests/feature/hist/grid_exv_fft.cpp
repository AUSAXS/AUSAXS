#if defined(POCKETFFT_AVAILABLE)

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/GridExvFFT.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/BinEstimate.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <grid/Grid.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <form_factor/FormFactorType.h>
#include <settings/All.h>

#include <iostream>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::hist;

namespace {
    hist::detail::CompactCoordinatesFF<false> as_coordinates(const std::vector<Vector3<double>>& points) {
        std::vector<AtomFF> atoms(points.size());
        std::ranges::transform(points, atoms.begin(), [] (const Vector3<double>& p) {
            return AtomFF{p, form_factor::form_factor_t::EXCLUDED_VOLUME};
        });
        return hist::detail::factory::construct_ff<false>(atoms);
    }

    // the pair loop the lattice transform replaces, counting ordered pairs and leaving out self-pairs
    WeightedDistribution1D pair_loop(
        const hist::detail::CompactCoordinatesFF<false>& data_i,
        const hist::detail::CompactCoordinatesFF<false>& data_j,
        bool same_set, unsigned int bin_count)
    {
        WeightedDistribution1D p(bin_count);
        for (int i = 0; i < data_i.size(); ++i) {
            for (int j = same_set ? i+1 : 0; j < data_j.size(); ++j) {
                hist::detail::evaluate1<false, 2>(p, data_i, data_j, i, j);
            }
        }
        return p;
    }

    void compare(const WeightedDistribution1D& expected, const WeightedDistribution1D& actual) {
        REQUIRE(expected.size() == actual.size());
        for (unsigned int i = 0; i < expected.size(); ++i) {
            if (expected.index(i).count != actual.index(i).count) {
                std::cout << "grid_exv_fft: bin " << i << " holds " << actual.index(i).count
                          << " pairs, expected " << expected.index(i).count << std::endl;
            }
            REQUIRE(expected.index(i).count == actual.index(i).count);

            // the transform forms its distances in double while the pair loop uses float, so they only agree to float precision
            REQUIRE_THAT(actual.index(i).bin_center, Catch::Matchers::WithinRel(expected.index(i).bin_center, 1e-6));
        }
    }

    double total_count(const WeightedDistribution1D& p) {
        double sum = 0;
        for (unsigned int i = 0; i < p.size(); ++i) {sum += p.index(i).count;}
        return sum;
    }
}

// The transform must reproduce the pair loop it replaces exactly - the counts are integers, so there is no tolerance
// to hide behind.
TEST_CASE("lattice::self_correlation: matches the pair loop", "[files]") {
    settings::general::verbose = false;
    settings::grid::cell_width = 1;
    settings::grid::exv::width = 1;

    Molecule protein("tests/files/2epe.pdb");
    protein.clear_hydration();
    auto exv = grid::exv::RawGridExv::create(protein.get_grid());
    REQUIRE(exv.spacing == 1);
    REQUIRE(!exv.interior.empty());

    auto data_x = as_coordinates(exv.interior);
    unsigned int bin_count = hist::detail::required_bin_count<false>(data_x);
    double inv_bin_width = hist::detail::WidthController<false>::get_inv_width();

    REQUIRE(exv.interior_sites.size() == exv.interior.size());

    auto lattice = hist::detail::lattice::self_correlation(exv, inv_bin_width, bin_count);
    auto n = static_cast<double>(exv.interior.size());
    CHECK(total_count(lattice) == n*(n-1));
    compare(pair_loop(data_x, data_x, true, bin_count), lattice);
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

    auto data_x_i = as_coordinates(exv.interior);
    auto data_x_s = as_coordinates(exv.surface);
    unsigned int bin_count = hist::detail::required_bin_count<false>(data_x_i, data_x_s);
    double inv_bin_width = hist::detail::WidthController<false>::get_inv_width();

    REQUIRE(exv.interior_sites.size() == exv.interior.size());
    REQUIRE(exv.surface_sites.size() == exv.surface.size());

    auto lattice = hist::detail::lattice::correlations(exv, inv_bin_width, bin_count);
    SECTION("interior") {compare(pair_loop(data_x_i, data_x_i, true, bin_count), lattice.first);}
    SECTION("surface")  {compare(pair_loop(data_x_s, data_x_s, true, bin_count), lattice.second);}
    SECTION("cross")    {compare(pair_loop(data_x_i, data_x_s, false, bin_count), lattice.cross);}
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
