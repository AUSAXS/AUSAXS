#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/PlanarDihedralSymmetry.h>
#include <data/symmetry/DihedralSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <math/MatrixUtils.h>

#include <cmath>
#include <memory>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::symmetry;

namespace {
    // positions of a single body point p under {original + all copies}, for a body centred at cm
    std::vector<Vector3<double>> copy_positions(const IPolyhedralSymmetry& s, Vector3<double> cm, Vector3<double> p) {
        std::vector<Vector3<double>> out = {p};
        for (int rep = 1; rep <= static_cast<int>(s.repetitions()); ++rep) {out.push_back(s.get_transform(cm, rep)(p));}
        return out;
    }

    // centres of mass of the 2n subunits: the image of the body centre cm under each copy. It is
    // these (not arbitrary atoms, which rotate out of plane) that a planar dihedral keeps coplanar.
    std::vector<Vector3<double>> copy_centres(const IPolyhedralSymmetry& s, Vector3<double> cm) {
        return copy_positions(s, cm, cm);
    }
}

TEST_CASE("PlanarDihedralSymmetry: same group as the general dihedral") {
    int n = GENERATE(2, 3, 4, 6);
    PlanarDihedralSymmetry planar(n);
    DihedralSymmetry general(n);
    CHECK(planar.repetitions() == general.repetitions());                       // same 2n copies
    CHECK(planar.internal_pair_schedule().size() == general.internal_pair_schedule().size()); // same reuse classes
}

TEST_CASE("PlanarDihedralSymmetry: exposes only the two in-plane translation parameters") {
    PlanarDihedralSymmetry s(3);
    CHECK(s.span_translation().size() == 2);
    CHECK(s.span_rotation().size() == 3); // the frame orientation is still fully free
}

// The defining property: every copy must be coplanar, and this must hold for *any* frame
// orientation. The offset is interpreted in the group frame, so even after the whole group is
// rotated the copies stay in one plane -- the plane through the centre normal to the principal axis.
TEST_CASE("PlanarDihedralSymmetry: copies are coplanar under an arbitrary frame orientation") {
    int n = GENERATE(2, 3, 4, 5, 6);
    PlanarDihedralSymmetry s(n);

    // a deliberately generic frame + a generic in-plane offset
    s.rotation = {0.7, -1.3, 0.4};
    s.translation = {1.5, -0.8, 0.0}; // third component is pinned/ignored

    auto pts = copy_centres(s, {0.4, 0.2, 0.9});
    REQUIRE(pts.size() == 2*static_cast<std::size_t>(n));

    // fit a plane through three of the points and check every point lies on it
    Vector3<double> normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
    REQUIRE(normal.magnitude() > 1e-9);
    normal = normal / normal.magnitude();
    for (const auto& p : pts) {
        CHECK_THAT((p - pts[0]).dot(normal), Catch::Matchers::WithinAbs(0.0, 1e-9));
    }

    // and that plane's normal is the principal axis F*z_hat: pinning the axial offset is what makes
    // the copies coplanar, so the plane must be perpendicular to the principal axis
    Vector3<double> principal = matrix::rotation_matrix<double>(s.rotation)*Vector3<double>{0, 0, 1};
    CHECK_THAT(std::abs(normal.dot(principal)), Catch::Matchers::WithinAbs(1.0, 1e-9));
}

TEST_CASE("PlanarDihedralSymmetry: writing the pinned axial parameter has no effect") {
    PlanarDihedralSymmetry s(4);
    s.rotation = {0.2, 0.5, -0.9};
    s.translation = {1.0, 2.0, 0.0};
    auto before = copy_positions(s, {0.5, -0.2, 0.4}, {0.3, 0.6, 0.1});

    s.translation = {1.0, 2.0, 7.3}; // axial component should be ignored entirely
    auto after = copy_positions(s, {0.5, -0.2, 0.4}, {0.3, 0.6, 0.1});

    REQUIRE(before.size() == after.size());
    for (std::size_t i = 0; i < before.size(); ++i) {
        CHECK((before[i] - after[i]).magnitude() < 1e-12);
    }
}

TEST_CASE("PlanarDihedralSymmetry: general dihedral is genuinely non-planar (contrast)") {
    // sanity check that the constraint is doing something: a general D_n with an axial offset
    // produces two stacked rings, which are not coplanar
    DihedralSymmetry s(3);
    s.translation = {1.5, 0.0, 2.0}; // non-zero axial component -> stacked rings (principal axis = z here)
    auto pts = copy_centres(s, {0.4, 0.2, 0.9});

    Vector3<double> normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
    normal = normal / normal.magnitude();
    bool all_coplanar = true;
    for (const auto& p : pts) {if (std::abs((p - pts[0]).dot(normal)) > 1e-6) {all_coplanar = false; break;}}
    CHECK_FALSE(all_coplanar);
}

TEST_CASE("PlanarDihedralSymmetry: name parsing") {
    CHECK(get("dp3") == type::dp3);
    CHECK(dynamic_cast<PlanarDihedralSymmetry*>(create("dp3").get()) != nullptr);
    CHECK(dynamic_cast<PlanarDihedralSymmetry*>(create("dp2").get()) != nullptr);
    // the plain dihedral name stays the general (non-planar) group
    CHECK(dynamic_cast<PlanarDihedralSymmetry*>(create("d3").get()) == nullptr);
}
