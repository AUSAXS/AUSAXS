// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/selection/ParameterMask.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <math/Vector3.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::selection;

namespace {
    // A symmetry delta with every optimizable parameter set to a non-zero value, so masking is
    // observable as a transition to zero. Composite symmetries carry no spans of their own, so we
    // reach the leaves the same way ParameterMask::apply does.
    std::unique_ptr<symmetry::ISymmetry> nonzero_delta(std::string_view name, double value) {
        auto sym = symmetry::create(name);
        symmetry::for_each_leaf(*sym, [value](symmetry::ISymmetry& leaf) {
            for (auto& t : leaf.span_translation()) {t = value;}
            for (auto& r : leaf.span_rotation()) {r = value;}
        });
        return sym;
    }

    bool translation_all_zero(symmetry::ISymmetry& sym) {
        bool zero = true;
        symmetry::for_each_leaf(sym, [&zero](symmetry::ISymmetry& leaf) {
            auto s = leaf.span_translation();
            zero = zero && std::all_of(s.begin(), s.end(), [](double v) {return v == 0;});
        });
        return zero;
    }

    bool rotation_all_zero(symmetry::ISymmetry& sym) {
        bool zero = true;
        symmetry::for_each_leaf(sym, [&zero](symmetry::ISymmetry& leaf) {
            auto s = leaf.span_rotation();
            zero = zero && std::all_of(s.begin(), s.end(), [](double v) {return v == 0;});
        });
        return zero;
    }

    // Three distinct symmetries, one of them composite so the leaf recursion is exercised too.
    parameter::BodyTransformParametersRelative make_params() {
        parameter::BodyTransformParametersRelative params;
        params.translation = Vector3<double>{1, 2, 3};
        params.rotation = Vector3<double>{4, 5, 6};
        params.symmetry_pars.emplace();
        params.symmetry_pars->emplace_back(nonzero_delta("c3", 1));
        params.symmetry_pars->emplace_back(nonzero_delta("p2-c3", 2));
        params.symmetry_pars->emplace_back(nonzero_delta("c2", 3));
        return params;
    }
}

TEST_CASE("ParameterMask::apply untargeted") {
    // Regression safety: with target_symmetry unset, every symmetry must be treated exactly alike.
    SECTION("all() leaves everything untouched") {
        auto params = make_params();
        ParameterMask::all().apply(params);

        CHECK(params.translation.has_value());
        CHECK(params.rotation.has_value());
        REQUIRE(params.symmetry_pars.has_value());
        for (auto& sym : params.symmetry_pars.value()) {
            CHECK_FALSE(translation_all_zero(*sym));
            CHECK_FALSE(rotation_all_zero(*sym));
        }
    }

    SECTION("real_only() drops the symmetry parameters entirely") {
        auto params = make_params();
        ParameterMask::real_only().apply(params);

        CHECK(params.translation.has_value());
        CHECK(params.rotation.has_value());
        CHECK_FALSE(params.symmetry_pars.has_value());
    }

    SECTION("symmetry_only() clears the pose but keeps every symmetry active") {
        auto params = make_params();
        ParameterMask::symmetry_only().apply(params);

        CHECK_FALSE(params.translation.has_value());
        CHECK_FALSE(params.rotation.has_value());
        REQUIRE(params.symmetry_pars.has_value());
        for (auto& sym : params.symmetry_pars.value()) {
            CHECK_FALSE(translation_all_zero(*sym));
            CHECK_FALSE(rotation_all_zero(*sym));
        }
    }

    SECTION("symmetry_only_trans() zeroes the rotation of every symmetry") {
        auto params = make_params();
        ParameterMask::symmetry_only_trans().apply(params);

        REQUIRE(params.symmetry_pars.has_value());
        for (auto& sym : params.symmetry_pars.value()) {
            CHECK_FALSE(translation_all_zero(*sym));
            CHECK(rotation_all_zero(*sym));
        }
    }

    SECTION("symmetry_only_axis() zeroes the translation of every symmetry") {
        auto params = make_params();
        ParameterMask::symmetry_only_axis().apply(params);

        REQUIRE(params.symmetry_pars.has_value());
        for (auto& sym : params.symmetry_pars.value()) {
            CHECK(translation_all_zero(*sym));
            CHECK_FALSE(rotation_all_zero(*sym));
        }
    }
}

namespace {
    // how a select strategy composes a targeted mask: the class of parameters comes from the mask strategy, the slot from the drawn target
    ParameterMask targeting(ParameterMask mask, unsigned int isymmetry) {
        mask.target_symmetry = isymmetry;
        return mask;
    }
}

TEST_CASE("ParameterMask::apply targeted") {
    SECTION("only the targeted symmetry survives") {
        for (unsigned int target : {0u, 1u, 2u}) {
            auto params = make_params();
            targeting(ParameterMask::symmetry_only(), target).apply(params);

            REQUIRE(params.symmetry_pars.has_value());
            REQUIRE(params.symmetry_pars->size() == 3);
            for (unsigned int i = 0; i < 3; ++i) {
                auto& sym = *params.symmetry_pars.value()[i];
                CHECK(translation_all_zero(sym) == (i != target));
                CHECK(rotation_all_zero(sym) == (i != target));
            }
        }
    }

    SECTION("targeting composes with a narrower symmetry mask") {
        // the slot and the parameter class are two independent decisions: the target says which symmetry moves, the mask which of its parameters may
        auto params = make_params();
        targeting(ParameterMask::symmetry_only_axis(), 1).apply(params);

        REQUIRE(params.symmetry_pars.has_value());
        REQUIRE(params.symmetry_pars->size() == 3);
        for (unsigned int i = 0; i < 3; ++i) {
            auto& sym = *params.symmetry_pars.value()[i];
            CHECK(translation_all_zero(sym));               // axis-only, so no symmetry translates
            CHECK(rotation_all_zero(sym) == (i != 1));      // and only the targeted one rotates
        }
    }

    SECTION("the host body's own pose is frozen") {
        auto params = make_params();
        targeting(ParameterMask::symmetry_only(), 1).apply(params);

        CHECK_FALSE(params.translation.has_value());
        CHECK_FALSE(params.rotation.has_value());
    }

    SECTION("targeting the only symmetry is equivalent to symmetry_only()") {
        parameter::BodyTransformParametersRelative targeted;
        targeted.symmetry_pars.emplace();
        targeted.symmetry_pars->emplace_back(nonzero_delta("c3", 1));
        targeting(ParameterMask::symmetry_only(), 0).apply(targeted);

        parameter::BodyTransformParametersRelative untargeted;
        untargeted.symmetry_pars.emplace();
        untargeted.symmetry_pars->emplace_back(nonzero_delta("c3", 1));
        ParameterMask::symmetry_only().apply(untargeted);

        REQUIRE(targeted.symmetry_pars.has_value());
        REQUIRE(untargeted.symmetry_pars.has_value());
        auto lhs_t = targeted.symmetry_pars.value()[0]->span_translation();
        auto rhs_t = untargeted.symmetry_pars.value()[0]->span_translation();
        auto lhs_r = targeted.symmetry_pars.value()[0]->span_rotation();
        auto rhs_r = untargeted.symmetry_pars.value()[0]->span_rotation();
        CHECK(std::equal(lhs_t.begin(), lhs_t.end(), rhs_t.begin(), rhs_t.end()));
        CHECK(std::equal(lhs_r.begin(), lhs_r.end(), rhs_r.begin(), rhs_r.end()));
    }
}
