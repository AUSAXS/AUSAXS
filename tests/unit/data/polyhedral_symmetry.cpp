#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/TetrahedralSymmetry.h>
#include <data/symmetry/OctahedralSymmetry.h>
#include <data/symmetry/IcosahedralSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include <functional>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    // a small, deliberately asymmetric test body
    const std::vector<Vector3<double>> body = {
        {1.0, 0.0, 0.0}, {0.3, 1.7, 0.2}, {-0.5, 0.4, 2.1}
    };
    const Vector3<double> cm = {0.0, 0.0, 0.0};

    // each concrete group, behind the shared PolyhedralSymmetry interface
    std::function<std::unique_ptr<IPolyhedralSymmetry>()> make_tetra = [] () {return std::make_unique<TetrahedralSymmetry>();};
    std::function<std::unique_ptr<IPolyhedralSymmetry>()> make_octa  = [] () {return std::make_unique<OctahedralSymmetry>();};
    std::function<std::unique_ptr<IPolyhedralSymmetry>()> make_icosa = [] () {return std::make_unique<IcosahedralSymmetry>();};

    // every distance between two placements of the body, brute-force
    std::vector<double> cross_distances(const IPolyhedralSymmetry& s, int repA, int repB) {
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

    // number of geometrically distinct points (the size of an orbit)
    int count_distinct(const std::vector<Vector3<double>>& pts) {
        int n = 0;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            bool seen = false;
            for (std::size_t j = 0; j < i; ++j) {if ((pts[i]-pts[j]).magnitude() < 1e-9) {seen = true; break;}}
            if (!seen) {++n;}
        }
        return n;
    }

    // number of distinct images of a single atom under {original + all copies}; equals
    // |group| / |stabilizer of p|, so it collapses when p lies on a symmetry axis
    int orbit_size(const IPolyhedralSymmetry& s, Vector3<double> p) {
        std::vector<Vector3<double>> orbit = {p}; // rep 0 = the original atom
        for (int rep = 1; rep <= static_cast<int>(s.repetitions()); ++rep) {orbit.push_back(s.get_transform({0, 0, 0}, rep)(p));}
        return count_distinct(orbit);
    }

    const double phi = (1 + std::sqrt(5.0))/2; // golden ratio; (0, 1, phi) is an icosahedral 5-fold axis
}

TEST_CASE("PolyhedralSymmetry: group order") {
    CHECK(TetrahedralSymmetry().repetitions() == 11);
    CHECK(OctahedralSymmetry().repetitions()  == 23);
    CHECK(IcosahedralSymmetry().repetitions() == 59);
}

TEST_CASE("PolyhedralSymmetry: pair schedule covers every copy-pair exactly once") {
    auto make = GENERATE(make_tetra, make_octa, make_icosa);
    auto s = make();
    int n = static_cast<int>(s->repetitions()) + 1; // bodies including the original

    long total = 0;
    for (const auto& pair : s->internal_pair_schedule()) {
        CHECK(0 <= pair.repA); CHECK(pair.repA < n);
        CHECK(0 <= pair.repB); CHECK(pair.repB < n);
        CHECK(0 < pair.scale);
        total += pair.scale;
    }
    CHECK(total == static_cast<long>(n)*(n-1)/2); // C(n,2): every unordered pair represented
}

TEST_CASE("PolyhedralSymmetry: schedule representatives reproduce all inter-copy distances") {
    auto make = GENERATE(make_tetra, make_octa);
    auto s = make();
    int n = static_cast<int>(s->repetitions()) + 1;

    // brute-force multiset of every inter-copy distance
    std::vector<double> brute;
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            auto d = cross_distances(*s, i, j);
            brute.insert(brute.end(), d.begin(), d.end());
        }
    }
    std::sort(brute.begin(), brute.end());

    // reconstruct it from the schedule: one representative per class, weighted by scale
    std::vector<double> reconstructed;
    for (const auto& pair : s->internal_pair_schedule()) {
        auto d = cross_distances(*s, pair.repA, pair.repB);
        for (int k = 0; k < pair.scale; ++k) {reconstructed.insert(reconstructed.end(), d.begin(), d.end());}
    }
    std::sort(reconstructed.begin(), reconstructed.end());

    REQUIRE(reconstructed.size() == brute.size());
    for (std::size_t k = 0; k < brute.size(); ++k) {
        CHECK_THAT(reconstructed[k], Catch::Matchers::WithinAbs(brute[k], 1e-6));
    }
}

TEST_CASE("PolyhedralSymmetry::get_transform: copies are proper rotations about the centre") {
    auto make = GENERATE(make_tetra, make_octa, make_icosa);
    auto s = make();

    SECTION("every copy is a proper rotation about the group centre") {
        for (int rep = 1; rep <= static_cast<int>(s->repetitions()); ++rep) {
            auto f = s->get_transform({0, 0, 0}, rep);
            // with no offset the group centre is the origin, which every rotation fixes
            CHECK(f({0, 0, 0}) == Vector3<double>(0, 0, 0));
            // a proper rotation preserves lengths
            CHECK_THAT(f({1, 0, 0}).magnitude(), Catch::Matchers::WithinAbs(1.0, 1e-9));
            CHECK_THAT(f({0.3, 1.7, 0.2}).magnitude(), Catch::Matchers::WithinAbs(Vector3<double>(0.3, 1.7, 0.2).magnitude(), 1e-9));
        }
    }

    SECTION("offset moves the fixed point to c = cm + translation") {
        s->translation = {1, 2, 3};
        auto f = s->get_transform({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 2, 3));
    }
}

// A single atom is replicated to one point per group element (|group| = 12 / 24 / 60), but when it
// sits on a symmetry axis the rotations about that axis map it to itself, collapsing the orbit. The
// surviving count is |group| / |stabilizer|, which is how the familiar low-vertex shapes appear.
const Vector3<double> generic{0.31, 0.59, 0.83}; // off every symmetry axis of every group

TEST_CASE("PolyhedralSymmetry: tetrahedral orbit sizes at special positions") {
    TetrahedralSymmetry s;
    CHECK(orbit_size(s, {1, 1, 1}) == 4);  // 3-fold axis -> 4 tetrahedron vertices
    CHECK(orbit_size(s, {0, 0, 1}) == 6);  // 2-fold axis -> 6 octahedron vertices
    CHECK(orbit_size(s, generic)   == 12); // generic position -> full 12-element orbit
}

TEST_CASE("PolyhedralSymmetry: octahedral orbit sizes at special positions") {
    OctahedralSymmetry s;
    CHECK(orbit_size(s, {0, 0, 1}) == 6);  // 4-fold axis -> 6 octahedron vertices
    CHECK(orbit_size(s, {1, 1, 1}) == 8);  // 3-fold axis -> 8 cube vertices
    CHECK(orbit_size(s, {1, 1, 0}) == 12); // 2-fold axis -> 12 edge midpoints (cuboctahedron)
    CHECK(orbit_size(s, generic)   == 24); // generic position -> full 24-element orbit
}

TEST_CASE("PolyhedralSymmetry: icosahedral orbit sizes at special positions") {
    IcosahedralSymmetry s;
    CHECK(orbit_size(s, {0, 1, phi}) == 12); // 5-fold axis -> 12 icosahedron vertices
    CHECK(orbit_size(s, {1, 1, 1})   == 20); // 3-fold axis -> 20 dodecahedron vertices
    CHECK(orbit_size(s, {1, 0, 0})   == 30); // 2-fold axis -> 30 edge midpoints
    CHECK(orbit_size(s, generic)     == 60); // generic position -> full 60-element orbit
}

TEST_CASE("PolyhedralSymmetry: a rigid line shows the rotation of each copy") {
    // Three collinear, equally-spaced atoms form a small "arrow" whose orientation is visible (unlike a
    // single point). Placed off every symmetry axis its orbit is the full group, which makes the rotation
    // of each copy explicit: every copy is the same arrow rigidly rotated to point a new way.
    auto make = GENERATE(make_tetra, make_octa, make_icosa);
    auto s = make();

    const Vector3<double> step{1.0, 0.5, 0.25};
    const Vector3<double> p0{0.3, 0.1, 0.7};
    const std::array<Vector3<double>, 3> line = {p0, p0 + step, p0 + 2*step};

    std::vector<Vector3<double>> heads = {line[0]}; // rep 0 = the original arrow
    for (int rep = 1; rep <= static_cast<int>(s->repetitions()); ++rep) {
        auto t = s->get_transform({0, 0, 0}, rep);
        Vector3<double> q0 = t(line[0]), q1 = t(line[1]), q2 = t(line[2]);

        // the copy is still a straight, equally-spaced line: the two steps stay identical
        // vectors, which only holds if the copy is a rigid rotation (no shear, no scaling)
        CHECK(((q1 - q0) - (q2 - q1)).magnitude() < 1e-9);
        // ... the step length is preserved
        CHECK_THAT((q1 - q0).magnitude(), Catch::Matchers::WithinAbs(step.magnitude(), 1e-9));
        // ... but the arrow now points a different way: the rotation genuinely reorients it
        CHECK((q1 - q0 - step).magnitude() > 1e-6);
        heads.emplace_back(q0);
    }

    // every rotation carries the arrow to a distinct position (full orbit of a generic body): 12 / 24 / 60
    CHECK(count_distinct(heads) == static_cast<int>(s->repetitions()) + 1);
}

TEST_CASE("PolyhedralSymmetry: name parsing maps to the right concrete group") {
    CHECK(dynamic_cast<TetrahedralSymmetry*>(create("t").get())           != nullptr);
    CHECK(dynamic_cast<TetrahedralSymmetry*>(create("tetrahedral").get()) != nullptr);
    CHECK(dynamic_cast<OctahedralSymmetry*>(create("o").get())            != nullptr);
    CHECK(dynamic_cast<OctahedralSymmetry*>(create("octahedral").get())   != nullptr);
    CHECK(dynamic_cast<IcosahedralSymmetry*>(create("i").get())           != nullptr);
    CHECK(dynamic_cast<IcosahedralSymmetry*>(create("icosahedral").get()) != nullptr);

    CHECK(get("t") == type::tetrahedral);
    CHECK(get("o") == type::octahedral);
    CHECK(get("i") == type::icosahedral);
    CHECK_THROWS(get("not-a-symmetry-name"));
}