#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::data;

TEST_CASE("UniformParameterGenerator::next") {
    settings::general::verbose = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    int iterations = 100;
    double length_start = GENERATE(1, 2, 3);
    double rad_start = GENERATE(1, 2, 3);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)})}});

    SECTION("RotationsOnly") {
        rigidbody::parameter::UniformParameterGenerator ro(&rb, iterations, {.rotation = rad_start});

        for (int i = 0; i < iterations; i++) {
            auto p = ro.next(0);
            CHECK_FALSE(p.translation.has_value());
            REQUIRE(-rad_start        <= p.rotation.value().x());
            REQUIRE(p.rotation.value().x()    <= rad_start     );
            REQUIRE(-rad_start        <= p.rotation.value().y());
            REQUIRE(p.rotation.value().y()    <= rad_start     );
            REQUIRE(-rad_start        <= p.rotation.value().z());
            REQUIRE(p.rotation.value().z()    <= rad_start     );
        }
    }

    SECTION("TranslationsOnly") {
        rigidbody::parameter::UniformParameterGenerator to(&rb, iterations, {.translation = length_start});

        for (int i = 0; i < iterations; i++) {
            auto p = to.next(0);
            REQUIRE(-length_start     <= p.translation.value().x());
            REQUIRE(p.translation.value().x() <= length_start     );
            REQUIRE(-length_start     <= p.translation.value().y());
            REQUIRE(p.translation.value().y() <= length_start     );
            REQUIRE(-length_start     <= p.translation.value().z());
            REQUIRE(p.translation.value().z() <= length_start     );
            CHECK_FALSE(p.rotation.has_value());
        }
    }

    SECTION("translations & rotations") {
        rigidbody::parameter::UniformParameterGenerator ap(&rb, iterations, {.translation = length_start, .rotation = rad_start});

        for (int i = 0; i < iterations; i++) {
            auto p = ap.next(0);
            REQUIRE(-length_start     <= p.translation.value().x());
            REQUIRE(p.translation.value().x() <= length_start     );
            REQUIRE(-length_start     <= p.translation.value().y());
            REQUIRE(p.translation.value().y() <= length_start     );
            REQUIRE(-length_start     <= p.translation.value().z());
            REQUIRE(p.translation.value().z() <= length_start     );
            REQUIRE(-rad_start        <= p.rotation.value().x()   );
            REQUIRE(p.rotation.value().x()    <= rad_start        );
            REQUIRE(-rad_start        <= p.rotation.value().y()   );
            REQUIRE(p.rotation.value().y()    <= rad_start        );
            REQUIRE(-rad_start        <= p.rotation.value().z()   );
            REQUIRE(p.rotation.value().z()    <= rad_start        );
        }
    }
}