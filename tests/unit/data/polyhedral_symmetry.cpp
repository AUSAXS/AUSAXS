#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/PolyhedralSymmetry.h>

#include <algorithm>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    // a small, deliberately asymmetric test body
    const std::vector<Vector3<double>> body = {
        {1.0, 0.0, 0.0}, {0.3, 1.7, 0.2}, {-0.5, 0.4, 2.1}
    };
    const Vector3<double> cm = {0.0, 0.0, 0.0};

    // every distance between two placements of the body, brute-force
    std::vector<double> cross_distances(const PolyhedralSymmetry& s, int repA, int repB) {
        auto build = [&](int rep) {
            std::vector<Vector3<double>> out;
            if (rep == 0) {out = body;}
            else {
                auto t = s.get_transform(cm, rep);
                for (const auto& v : body) {out.push_back(t(v));}
            }
            return out;
        };
        auto A = build(repA), B = build(repB);
        std::vector<double> d;
        for (const auto& a : A) {for (const auto& b : B) {d.push_back((a-b).magnitude());}}
        std::sort(d.begin(), d.end());
        return d;
    }
}

TEST_CASE("PolyhedralSymmetry: group order") {
    CHECK(PolyhedralSymmetry(PolyhedralGroup::tetrahedral).repetitions() == 11);
    CHECK(PolyhedralSymmetry(PolyhedralGroup::octahedral).repetitions()  == 23);
    CHECK(PolyhedralSymmetry(PolyhedralGroup::icosahedral).repetitions() == 59);
}

TEST_CASE("PolyhedralSymmetry: pair schedule covers every copy-pair exactly once") {
    auto group = GENERATE(PolyhedralGroup::tetrahedral, PolyhedralGroup::octahedral, PolyhedralGroup::icosahedral);
    PolyhedralSymmetry s(group);
    int n = static_cast<int>(s.repetitions()) + 1; // bodies including the original

    long total = 0;
    for (const auto& pair : s.internal_pair_schedule()) {
        CHECK(0 <= pair.repA); CHECK(pair.repA < n);
        CHECK(0 <= pair.repB); CHECK(pair.repB < n);
        CHECK(0 < pair.scale);
        total += pair.scale;
    }
    CHECK(total == static_cast<long>(n)*(n-1)/2); // C(n,2): every unordered pair represented
}

TEST_CASE("PolyhedralSymmetry: schedule representatives reproduce all inter-copy distances") {
    auto group = GENERATE(PolyhedralGroup::tetrahedral, PolyhedralGroup::octahedral);
    PolyhedralSymmetry s(group);
    int n = static_cast<int>(s.repetitions()) + 1;

    // brute-force multiset of every inter-copy distance
    std::vector<double> brute;
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            auto d = cross_distances(s, i, j);
            brute.insert(brute.end(), d.begin(), d.end());
        }
    }
    std::sort(brute.begin(), brute.end());

    // reconstruct it from the schedule: one representative per class, weighted by scale
    std::vector<double> reconstructed;
    for (const auto& pair : s.internal_pair_schedule()) {
        auto d = cross_distances(s, pair.repA, pair.repB);
        for (int k = 0; k < pair.scale; ++k) {reconstructed.insert(reconstructed.end(), d.begin(), d.end());}
    }
    std::sort(reconstructed.begin(), reconstructed.end());

    REQUIRE(reconstructed.size() == brute.size());
    for (std::size_t k = 0; k < brute.size(); ++k) {
        CHECK_THAT(reconstructed[k], Catch::Matchers::WithinAbs(brute[k], 1e-6));
    }
}
