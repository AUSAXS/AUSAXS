#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>

#include <algorithm>
#include <numbers>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    std::unique_ptr<CyclicSymmetry> cyclic(double angle, int reps, Vector3<double> offset) {
        return std::make_unique<CyclicSymmetry>(
            CyclicSymmetry::_Relation{offset}, CyclicSymmetry::_Repeat{{0, 0, 1}, angle}, reps
        );
    }
}

TEST_CASE("CompositeSymmetry: repetition count") {
    // p2 (1 copy) nested in c3 (2 copies) -> 2*3 = 6 placements -> 5 repetitions
    CompositeSymmetry p2_c3(
        std::make_unique<PointSymmetry>(Vector3<double>{3, 0, 0}, Vector3<double>{0, 0, 0}),
        cyclic(2*std::numbers::pi/3, 2, {5, 0, 0})
    );
    CHECK(p2_c3.repetitions() == 5);

    // c3 (2) nested in c3 (2) -> 3*3 = 9 placements -> 8 repetitions
    CompositeSymmetry c3_c3(
        cyclic(2*std::numbers::pi/3, 2, {3, 0, 0}),
        cyclic(2*std::numbers::pi/3, 2, {7, 0, 0})
    );
    CHECK(c3_c3.repetitions() == 8);
}

TEST_CASE("CompositeSymmetry: pair schedule covers every copy-pair exactly once") {
    CompositeSymmetry sym(
        std::make_unique<PointSymmetry>(Vector3<double>{4, 1, 0}, Vector3<double>{0, 0, 0}),
        cyclic(2*std::numbers::pi/3, 2, {6, 0, 0})
    );
    int n = static_cast<int>(sym.repetitions()) + 1;

    long total = 0;
    for (const auto& pair : sym.internal_pair_schedule()) {
        CHECK(0 <= pair.repA); CHECK(pair.repA < n);
        CHECK(0 <= pair.repB); CHECK(pair.repB < n);
        CHECK(0 < pair.scale);
        total += pair.scale;
    }
    CHECK(total == static_cast<long>(n)*(n-1)/2);
}

TEST_CASE("CompositeSymmetry: schedule reproduces all inter-copy distances") {
    const std::vector<Vector3<double>> body = {{1.0, 0.0, 0.0}, {0.3, 1.7, 0.2}, {-0.5, 0.4, 2.1}};
    const Vector3<double> cm = {0.0, 0.0, 0.0};

    CompositeSymmetry sym(
        cyclic(std::numbers::pi, 1, {3, 0, 0}),          // inner c2
        cyclic(2*std::numbers::pi/3, 2, {7, 0, 0})       // outer c3
    );
    int n = static_cast<int>(sym.repetitions()) + 1;

    auto placement = [&](int rep) {
        std::vector<Vector3<double>> out;
        if (rep == 0) {return body;}
        auto t = sym.get_transform(cm, rep);
        for (const auto& v : body) {out.push_back(t(v));}
        return out;
    };
    auto cross = [&](int a, int b) {
        auto A = placement(a), B = placement(b);
        std::vector<double> d;
        for (const auto& x : A) {for (const auto& y : B) {d.push_back((x-y).magnitude());}}
        return d;
    };

    std::vector<double> brute, reconstructed;
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            auto d = cross(i, j);
            brute.insert(brute.end(), d.begin(), d.end());
        }
    }
    for (const auto& pair : sym.internal_pair_schedule()) {
        auto d = cross(pair.repA, pair.repB);
        for (int k = 0; k < pair.scale; ++k) {reconstructed.insert(reconstructed.end(), d.begin(), d.end());}
    }
    std::sort(brute.begin(), brute.end());
    std::sort(reconstructed.begin(), reconstructed.end());

    REQUIRE(reconstructed.size() == brute.size());
    for (std::size_t k = 0; k < brute.size(); ++k) {
        CHECK_THAT(reconstructed[k], Catch::Matchers::WithinAbs(brute[k], 1e-6));
    }
}
