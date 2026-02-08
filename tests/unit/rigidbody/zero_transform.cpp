#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("ZeroTransform: zeroed delta parameters do not cause drift") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a3({0, 5, 0}, form_factor::form_factor_t::C);
    Body b1(std::vector<AtomFF>{a1});
    Body b2(std::vector<AtomFF>{a2});
    Body b3(std::vector<AtomFF>{a3});
    Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3}});

    auto& transformer = rigidbody.transformer;

    SECTION("repeated zero transforms do not drift rotation") {
        auto original_rotation = rigidbody.conformation->absolute_parameters.parameters[0].rotation;
        auto original_translation = rigidbody.conformation->absolute_parameters.parameters[0].translation;

        for (int i = 0; i < 1000; ++i) {
            parameter::BodyTransformParametersRelative zero_params({0, 0, 0}, {0, 0, 0});
            transformer->apply(std::move(zero_params), 0u);
        }

        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(params.rotation.x(), Catch::Matchers::WithinAbs(original_rotation.x(), 1e-10));
        REQUIRE_THAT(params.rotation.y(), Catch::Matchers::WithinAbs(original_rotation.y(), 1e-10));
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(original_rotation.z(), 1e-10));
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(original_translation.x(), 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(original_translation.y(), 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(original_translation.z(), 1e-10));
    }

    SECTION("zero transforms after non-zero transform preserve state") {
        // First apply a non-trivial transform
        parameter::BodyTransformParametersRelative real_params({1.5, 2.3, -0.7}, {0.3, -0.1, 0.5});
        transformer->apply(std::move(real_params), 0u);

        auto rotation_after = rigidbody.conformation->absolute_parameters.parameters[0].rotation;
        auto translation_after = rigidbody.conformation->absolute_parameters.parameters[0].translation;

        // Then apply many zero transforms
        for (int i = 0; i < 1000; ++i) {
            parameter::BodyTransformParametersRelative zero_params({0, 0, 0}, {0, 0, 0});
            transformer->apply(std::move(zero_params), 0u);
        }

        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(params.rotation.x(), Catch::Matchers::WithinAbs(rotation_after.x(), 1e-10));
        REQUIRE_THAT(params.rotation.y(), Catch::Matchers::WithinAbs(rotation_after.y(), 1e-10));
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(rotation_after.z(), 1e-10));
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(translation_after.x(), 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(translation_after.y(), 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(translation_after.z(), 1e-10));
    }
}
