#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/DihedralSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    // a small, deliberately asymmetric test body
    const std::vector<Vector3<double>> body = {
        {1.0, 0.0, 0.0}, {0.3, 1.7, 0.2}, {-0.5, 0.4, 2.1}
    };
    const Vector3<double> cm = {0.0, 0.0, 0.0};

    std::vector<double> cross_distances(const IPolyhedralSymmetry& s, int repA, int repB) {
        auto build = [&](int rep) {
            std::vector<Vector3<double>> out;
            if (rep == 0) {out = body;}
            else {auto t = s.get_transform(cm, rep); for (const auto& v : body) {out.push_back(t(v));}}
            return out;
        };
        auto A = build(repA), B = build(repB);
        std::vector<double> d;
        for (const auto& a : A) {for (const auto& b : B) {d.push_back((a-b).magnitude());}}
        std::sort(d.begin(), d.end());
        return d;
    }

    int count_distinct(const std::vector<Vector3<double>>& pts) {
        int n = 0;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            bool seen = false;
            for (std::size_t j = 0; j < i; ++j) {if ((pts[i]-pts[j]).magnitude() < 1e-9) {seen = true; break;}}
            if (!seen) {++n;}
        }
        return n;
    }

    int orbit_size(const IPolyhedralSymmetry& s, Vector3<double> p) {
        std::vector<Vector3<double>> orbit = {p};
        for (int rep = 1; rep <= static_cast<int>(s.repetitions()); ++rep) {orbit.push_back(s.get_transform({0, 0, 0}, rep)(p));}
        return count_distinct(orbit);
    }
}

TEST_CASE("DihedralSymmetry: group order is 2n") {
    for (int n = 2; n <= 12; ++n) {
        DihedralSymmetry s(n);
        CHECK(s.repetitions() == static_cast<unsigned int>(2*n) - 1); // 2n copies including the original
    }
}

TEST_CASE("DihedralSymmetry: pair schedule covers every copy-pair exactly once") {
    int n = GENERATE(2, 3, 4, 5, 6, 8, 12);
    DihedralSymmetry s(n);
    int m = static_cast<int>(s.repetitions()) + 1;

    long total = 0;
    for (const auto& pair : s.internal_pair_schedule()) {
        CHECK(0 <= pair.repA); CHECK(pair.repA < m);
        CHECK(0 <= pair.repB); CHECK(pair.repB < m);
        CHECK(0 < pair.scale);
        total += pair.scale;
    }
    CHECK(total == static_cast<long>(m)*(m-1)/2); // C(m,2): every unordered pair represented once
}

TEST_CASE("DihedralSymmetry: schedule representatives reproduce all inter-copy distances") {
    // the whole reason the group is modelled explicitly: one representative per equal-distance class,
    // weighted by scale, must reconstruct the full inter-copy distance multiset
    int n = GENERATE(2, 3, 4, 6);
    DihedralSymmetry s(n);
    int m = static_cast<int>(s.repetitions()) + 1;

    std::vector<double> brute;
    for (int i = 0; i < m; ++i) {
        for (int j = i+1; j < m; ++j) {
            auto d = cross_distances(s, i, j);
            brute.insert(brute.end(), d.begin(), d.end());
        }
    }
    std::sort(brute.begin(), brute.end());

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

TEST_CASE("DihedralSymmetry: copies are proper rotations about the centre") {
    int n = GENERATE(2, 3, 5, 6);
    DihedralSymmetry s(n);
    for (int rep = 1; rep <= static_cast<int>(s.repetitions()); ++rep) {
        auto f = s.get_transform({0, 0, 0}, rep);
        CHECK(f({0, 0, 0}) == Vector3<double>(0, 0, 0));               // centre is fixed
        CHECK_THAT(f({1, 0, 0}).magnitude(), Catch::Matchers::WithinAbs(1.0, 1e-9)); // length preserved
    }
}

// The defining property: D_n has n two-fold axes *perpendicular* to the principal axis, so it is not
// the cyclic C_2n. On the principal axis the C_n rotations act trivially but the perpendicular flips
// send the axis to its negative, giving an orbit of exactly two points (a single fixed point would
// mean the group was cyclic about that axis).
TEST_CASE("DihedralSymmetry: perpendicular two-fold axes distinguish D_n from C_2n") {
    int n = GENERATE(2, 3, 4, 6);
    DihedralSymmetry s(n);

    SECTION("principal axis maps to its negative under some copy (not fixed as in C_2n)") {
        CHECK(orbit_size(s, {0, 0, 1}) == 2);
        bool found_flip = false;
        for (int rep = 1; rep <= static_cast<int>(s.repetitions()); ++rep) {
            if ((s.get_transform({0, 0, 0}, rep)({0, 0, 1}) - Vector3<double>{0, 0, -1}).magnitude() < 1e-9) {found_flip = true;}
        }
        CHECK(found_flip);
    }

    SECTION("a generic point has the full 2n-element orbit") {
        CHECK(orbit_size(s, {0.31, 0.59, 0.83}) == 2*n);
    }
}

TEST_CASE("DihedralSymmetry: name parsing") {
    CHECK(get("d3") == type::d3);
    CHECK(dynamic_cast<DihedralSymmetry*>(create("d3").get()) != nullptr);

    // a c2 + single c_n pair (either order) is the dihedral group, not a generic composite
    for (const auto& name : {"c2-c3", "c3-c2"}) {
        auto s = create(name);
        auto d = dynamic_cast<DihedralSymmetry*>(s.get());
        REQUIRE(d != nullptr);
        CHECK(d->repetitions() == 5); // D3: six copies
    }
    CHECK(dynamic_cast<DihedralSymmetry*>(create("c2-c2").get()) != nullptr); // D2

    // deeper chains and non-c2 pairs remain generic composites
    CHECK(dynamic_cast<CompositeSymmetry*>(create("c2-c2-c3").get()) != nullptr);
    CHECK(dynamic_cast<CompositeSymmetry*>(create("c3-c4").get())    != nullptr);
}
