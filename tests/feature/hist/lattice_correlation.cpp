#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/LatticeCorrelation.h>
#include <hist/detail/CompactCoordinatesFF.h>
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
        std::transform(points.begin(), points.end(), atoms.begin(), [] (const Vector3<double>& p) {
            return AtomFF{p, form_factor::form_factor_t::EXCLUDED_VOLUME};
        });
        return hist::detail::CompactCoordinatesFF<false>(std::move(atoms));
    }

    // the pair loop the lattice transform replaces, counting ordered pairs and leaving out self-pairs
    WeightedDistribution1D pair_loop(
        const hist::detail::CompactCoordinatesFF<false>& data_i,
        const hist::detail::CompactCoordinatesFF<false>& data_j,
        bool same_set, unsigned int bin_count)
    {
        WeightedDistribution1D p(bin_count);
        for (int i = 0; i < static_cast<int>(data_i.size()); ++i) {
            for (int j = same_set ? i+1 : 0; j < static_cast<int>(data_j.size()); ++j) {
                hist::detail::evaluate1<false, 2>(p, data_i, data_j, i, j);
            }
        }
        return p;
    }

    void compare(const WeightedDistribution1D& expected, const WeightedDistribution1D& actual) {
        REQUIRE(expected.size() == actual.size());
        for (unsigned int i = 0; i < expected.size(); ++i) {
            if (expected.index(i).count != actual.index(i).count) {
                std::cout << "lattice_correlation: bin " << i << " holds " << actual.index(i).count
                          << " pairs, expected " << expected.index(i).count << std::endl;
            }
            REQUIRE(expected.index(i).count == actual.index(i).count);

            // the transform works from exact integer displacements where the pair loop works in float32, so the
            // accumulated distances agree only to the precision of the latter
            REQUIRE_THAT(
                actual.index(i).bin_center,
                Catch::Matchers::WithinRel(expected.index(i).bin_center, 1e-6) || Catch::Matchers::WithinAbs(0, 1e-6)
            );
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

    auto lattice = hist::detail::lattice::self_correlation(exv.interior, exv.spacing, inv_bin_width, bin_count);
    REQUIRE(lattice.has_value());

    const double n = static_cast<double>(exv.interior.size());
    CHECK(total_count(*lattice) == n*(n-1));
    compare(pair_loop(data_x, data_x, true, bin_count), *lattice);
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

    auto lattice = hist::detail::lattice::correlations(exv.interior, exv.surface, exv.spacing, inv_bin_width, bin_count);
    REQUIRE(lattice.has_value());

    SECTION("interior") {compare(pair_loop(data_x_i, data_x_i, true, bin_count), lattice->first);}
    SECTION("surface")  {compare(pair_loop(data_x_s, data_x_s, true, bin_count), lattice->second);}
    SECTION("cross")    {compare(pair_loop(data_x_i, data_x_s, false, bin_count), lattice->cross);}
}

// Everything must report failure rather than allocate past its budget or make up an answer for points that are not
// lattice-supported, since the callers rely on that to fall back to the pair loop.
TEST_CASE("lattice: refuses what it cannot do") {
    settings::general::verbose = false;
    unsigned int budget = settings::grid::exv::max_transform_memory;

    std::vector<Vector3<double>> lattice_points;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {lattice_points.emplace_back(i, j, k);}
        }
    }

    SECTION("lattice-supported points are accepted") {
        CHECK(hist::detail::lattice::self_correlation(lattice_points, 1, 1, 100).has_value());
    }

    SECTION("an unknown spacing is rejected") {
        CHECK(!hist::detail::lattice::self_correlation(lattice_points, 0, 1, 100).has_value());
    }

    SECTION("an empty set is rejected") {
        CHECK(!hist::detail::lattice::self_correlation({}, 1, 1, 100).has_value());
    }

    SECTION("off-lattice points are rejected") {
        auto off_lattice = lattice_points;
        off_lattice.back().z() += 0.5;
        CHECK(!hist::detail::lattice::self_correlation(off_lattice, 1, 1, 100).has_value());
    }

    SECTION("a wrong spacing is rejected") {
        CHECK(!hist::detail::lattice::self_correlation(lattice_points, 1.5, 1, 100).has_value());
    }

    SECTION("exceeding the memory budget is rejected") {
        settings::grid::exv::max_transform_memory = 0;
        CHECK(!hist::detail::lattice::self_correlation(lattice_points, 1, 1, 100).has_value());
        settings::grid::exv::max_transform_memory = budget;
    }
}
