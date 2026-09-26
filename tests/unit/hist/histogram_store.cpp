#include <catch2/catch_test_macros.hpp>

#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/data/WidthControllers.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <settings/GeneralSettings.h>

#include <vector>

using namespace ausaxs;
using namespace ausaxs::hist;

namespace {
    using Coordinates = hist::detail::CompactCoordinates<false>;
    using Calculator = distance_calculator::Calculator<false, false, UNIT_WEIGHTS>;

    // points on a line, each at the centre of the given distance bin from the origin
    Coordinates points_at(const std::vector<int>& bins) {
        double width = 1./hist::detail::WidthController<false>::get_inv_width();
        std::vector<Vector3<double>> points;
        points.reserve(bins.size());
        for (int bin : bins) {points.emplace_back(bin*width, 0, 0);}
        return hist::detail::factory::construct<false>(points);
    }

    void self(distance_calculator::HistogramStore<false>& store, const Coordinates& a, int id) {
        Calculator calculator(store);
        calculator.enqueue_calculate_self(a, id);
        calculator.run();
    }
}

TEST_CASE("HistogramStore") {
    settings::general::gpu = false;
    distance_calculator::HistogramStore<false> store(20, 2);

    SECTION("results are zero before the first run") {
        int id = store.allocate_1d();
        for (int i = 0; i < store.bins(); ++i) {CHECK(store.get_1d(id).index(i) == 0);}
    }

    SECTION("a run replaces the previous contents of a result") {
        int id = store.allocate_1d();
        auto a = points_at({0, 5});
        for (int run = 0; run < 2; ++run) {
            self(store, a, id);
            CHECK(store.get_1d(id).index(0) == 2); // each point with itself
            CHECK(store.get_1d(id).index(5) == 2); // the pair, counted in either order
        }
    }

    SECTION("a result calculated from an empty set is zeroed") {
        int id = store.allocate_1d();
        self(store, points_at({0, 5}), id);
        self(store, Coordinates{}, id);
        for (int i = 0; i < store.bins(); ++i) {CHECK(store.get_1d(id).index(i) == 0);}
    }

    SECTION("jobs into the same result sum within a run") {
        int id = store.allocate_1d();
        auto origin = points_at({0}), near = points_at({3}), far = points_at({7}); // kept alive until the calculator has run
        Calculator calculator(store);
        calculator.enqueue_calculate_cross(origin, near, id, 1);
        calculator.enqueue_calculate_cross(origin, far, id, 2);
        calculator.run();
        CHECK(store.get_1d(id).index(3) == 1);
        CHECK(store.get_1d(id).index(7) == 2);
    }

    SECTION("partitioned sets fill one histogram per class, and an empty class is zero") {
        int id = store.allocate_2d();
        std::vector<Coordinates> partitioned{points_at({0}), Coordinates{}};
        auto flat = points_at({4});
        Calculator calculator(store);
        calculator.enqueue_calculate_cross(partitioned, flat, id, 1);
        calculator.run();
        const auto& result = store.get_2d(id);
        CHECK(result.index(0, 4) == 1);
        for (int i = 0; i < store.bins(); ++i) {CHECK(result.index(1, i) == 0);}
    }
}
