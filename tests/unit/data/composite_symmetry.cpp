#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/Body.h>

#include <algorithm>
#include <numbers>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::symmetry;
using namespace ausaxs::data;

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

TEST_CASE("CompositeSymmetry: parameters are reached via for_each_leaf, not span_*") {
    CompositeSymmetry sym(
        cyclic(std::numbers::pi, 1, {3, 0, 0}),    // inner c2
        cyclic(2*std::numbers::pi/3, 2, {7, 0, 0}) // outer c3
    );

    // a composite has no single contiguous parameter span; calling span_*() directly is an error
    CHECK_THROWS(sym.span_translation());
    CHECK_THROWS(sym.span_rotation());

    // for_each_leaf instead visits the two sub-symmetries
    std::vector<ISymmetry*> leaves;
    for_each_leaf(sym, [&](ISymmetry& leaf) {leaves.push_back(&leaf);});
    REQUIRE(leaves.size() == 2);
    CHECK(leaves[0] == sym.inner.get());
    CHECK(leaves[1] == sym.outer.get());
    // and each leaf does expose a usable span
    CHECK(leaves[0]->span_translation().size() == 3);
}

namespace {
    // every placement of a test point under {original + all copies}
    std::vector<Vector3<double>> placements(ISymmetry& s, Vector3<double> cm, Vector3<double> p) {
        std::vector<Vector3<double>> out = {p};
        for (int rep = 1; rep <= static_cast<int>(s.repetitions()); ++rep) {out.push_back(s.get_transform(cm, rep)(p));}
        return out;
    }

    bool same_set(std::vector<Vector3<double>> a, std::vector<Vector3<double>> b) {
        if (a.size() != b.size()) {return false;}
        for (const auto& x : a) {
            auto it = std::find_if(b.begin(), b.end(), [&](const auto& y){return (x-y).magnitude() < 1e-9;});
            if (it == b.end()) {return false;}
            b.erase(it);
        }
        return b.empty();
    }
}

// Independent check of get_transform's placements: a composite of coprime cyclic symmetries about a shared axis/centre is the cyclic group 
// of their product order (C_m x C_n = C_{mn} iff gcd(m,n)=1), so its placement set must coincide with a single, separately-trusted c_{mn}
TEST_CASE("CompositeSymmetry: coprime cyclics about a shared axis reproduce a higher cyclic group") {
    const Vector3<double> cm{0, 0, 0};
    const Vector3<double> p{1.3, -0.7, 2.0};
    constexpr double tau = 2*std::numbers::pi;

    SECTION("c2 x c3 == c6") {
        CompositeSymmetry comp(cyclic(tau/2, 1, {0, 0, 0}), cyclic(tau/3, 2, {0, 0, 0}));
        auto c6 = cyclic(tau/6, 5, {0, 0, 0});
        REQUIRE(comp.repetitions() == c6->repetitions()); // 5
        CHECK(same_set(placements(comp, cm, p), placements(*c6, cm, p)));
    }

    SECTION("c3 x c4 == c12") {
        CompositeSymmetry comp(cyclic(tau/3, 2, {0, 0, 0}), cyclic(tau/4, 3, {0, 0, 0}));
        auto c12 = cyclic(tau/12, 11, {0, 0, 0});
        REQUIRE(comp.repetitions() == c12->repetitions()); // 11
        CHECK(same_set(placements(comp, cm, p), placements(*c12, cm, p)));
    }
}

// The coprime check above is blind to composition order (coaxial rotations commute). This pins the order by composing the *trusted* sub-symmetry 
// transforms by hand with non-commuting parts (a translation and a rotation), then checking the composite reproduces every outer(inner(p)).
TEST_CASE("CompositeSymmetry: get_transform composes outer after inner") {
    const Vector3<double> cm{0, 0, 0};
    const Vector3<double> p{1.3, -0.7, 2.0};

    auto inner_ref = std::make_unique<PointSymmetry>(Vector3<double>{10, 0, 0}, Vector3<double>{0, 0, 0}); // translation
    auto outer_ref = cyclic(std::numbers::pi, 1, {0, 0, 0});                                               // c2 rotation
    CompositeSymmetry comp(inner_ref->clone(), outer_ref->clone());

    auto apply = [&](ISymmetry& s, int idx, Vector3<double> v) {return idx == 0 ? v : s.get_transform(cm, idx)(v);};
    std::vector<Vector3<double>> reference; // outer(k) after inner(j) over the whole grid, identity at index 0
    for (int k = 0; k <= static_cast<int>(outer_ref->repetitions()); ++k) {
        for (int j = 0; j <= static_cast<int>(inner_ref->repetitions()); ++j) {
            reference.push_back(apply(*outer_ref, k, apply(*inner_ref, j, p)));
        }
    }
    CHECK(same_set(placements(comp, cm, p), reference));
}

// Nesting recursion: Composite(A, Composite(B, C)) must place every point at C(B(A(p))) over the full
// (A,B,C) grid. Built from non-commuting parts and checked against the *trusted* sub-transforms composed
// by hand, this independently validates the recursive decode/composition that the reused-vs-explicit
// feature test cannot (both of its paths share get_transform).
TEST_CASE("CompositeSymmetry: get_transform composes a 3-level nesting") {
    const Vector3<double> cm{0, 0, 0};
    const Vector3<double> p{1.3, -0.7, 2.0};

    auto A = std::make_unique<PointSymmetry>(Vector3<double>{10, 0, 0}, Vector3<double>{0, 0, 0}); // innermost: translation
    auto B = cyclic(std::numbers::pi, 1, {3, 0, 0});                                               // c2
    auto C = cyclic(2*std::numbers::pi/3, 2, {7, 0, 0});                                           // outermost: c3
    CompositeSymmetry comp(A->clone(), std::make_unique<CompositeSymmetry>(B->clone(), C->clone()));
    REQUIRE(comp.repetitions() == 11); // (1+1)(1+1)(1+2) - 1 = 11

    auto apply = [&](ISymmetry& s, int idx, Vector3<double> v) {return idx == 0 ? v : s.get_transform(cm, idx)(v);};
    std::vector<Vector3<double>> reference; // C(B(A(p))) over every (a, b, c) copy, identity at index 0
    for (int kc = 0; kc <= static_cast<int>(C->repetitions()); ++kc) {
        for (int jb = 0; jb <= static_cast<int>(B->repetitions()); ++jb) {
            for (int ia = 0; ia <= static_cast<int>(A->repetitions()); ++ia) {
                reference.push_back(apply(*C, kc, apply(*B, jb, apply(*A, ia, p))));
            }
        }
    }
    CHECK(same_set(placements(comp, cm, p), reference));
}

// A small concrete atomic system: explicit_structure must lay out na*(1+repetitions) atoms as [original, copy_1, ...], each block being
// get_transform(cm, rep) applied to the base atoms.
TEST_CASE("CompositeSymmetry: explicit_structure materialises every copy") {
    std::vector<AtomFF> base = {
        AtomFF({1, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({0, 2, 0}, form_factor::form_factor_t::C)
    };
    Body body{base};

    auto comp = std::make_unique<CompositeSymmetry>(
        std::make_unique<PointSymmetry>(Vector3<double>{10, 0, 0}, Vector3<double>{0, 0, 0}), // inner: translation
        cyclic(std::numbers::pi, 1, {0, 0, 0})                                                // outer: c2
    );
    auto* sym = comp.get();
    auto cm = body.get_cm();
    body.symmetry().add(std::move(comp));

    auto s = body.symmetry().explicit_structure();
    int na = static_cast<int>(base.size());
    REQUIRE(static_cast<int>(s.atoms.size()) == (1 + static_cast<int>(sym->repetitions()))*na); // 4 placements * 2 atoms

    for (int rep = 0; rep <= static_cast<int>(sym->repetitions()); ++rep) {
        for (int i = 0; i < na; ++i) {
            Vector3<double> expected = rep == 0 ? base[i].coordinates() : sym->get_transform(cm, rep)(base[i].coordinates());
            CHECK((s.atoms[rep*na + i].coordinates() - expected).magnitude() < 1e-9);
        }
    }
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
