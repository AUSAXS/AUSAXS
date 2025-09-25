#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/ParameterGenerationStrategies.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::data;

TEST_CASE("RotationsOnly::next") {
    settings::general::verbose = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    int iterations = 100;
    double length_start = GENERATE(1, 2, 3);
    double rad_start = GENERATE(1, 2, 3);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)})}});

    SECTION("RotationsOnly") {
        rigidbody::parameter::RotationsOnly ro(&rb, iterations, length_start, rad_start);

        for (int i = 0; i < iterations; i++) {
            auto p = ro.next(0);
            REQUIRE(p.translation.x() == 0             );
            REQUIRE(p.translation.y() == 0             );
            REQUIRE(p.translation.z() == 0             );
            REQUIRE(-rad_start        <= p.rotation.x());
            REQUIRE(p.rotation.x()    <= rad_start     );
            REQUIRE(-rad_start        <= p.rotation.y());
            REQUIRE(p.rotation.y()    <= rad_start     );
            REQUIRE(-rad_start        <= p.rotation.z());
            REQUIRE(p.rotation.z()    <= rad_start     );
        }
    }

    SECTION("TranslationsOnly") {
        rigidbody::parameter::TranslationsOnly to(&rb, iterations, length_start, rad_start);

        for (int i = 0; i < iterations; i++) {
            auto p = to.next(0);
            REQUIRE(-length_start     <= p.translation.x());
            REQUIRE(p.translation.x() <= length_start     );
            REQUIRE(-length_start     <= p.translation.y());
            REQUIRE(p.translation.y() <= length_start     );
            REQUIRE(-length_start     <= p.translation.z());
            REQUIRE(p.translation.z() <= length_start     );
            REQUIRE(p.rotation.x()    == 0                );
            REQUIRE(p.rotation.y()    == 0                );
            REQUIRE(p.rotation.z()    == 0                );
        }
    }

    SECTION("translations & rotations") {
        rigidbody::parameter::AllParameters ap(&rb, iterations, length_start, rad_start);

        for (int i = 0; i < iterations; i++) {
            auto p = ap.next(0);
            REQUIRE(-length_start     <= p.translation.x());
            REQUIRE(p.translation.x() <= length_start     );
            REQUIRE(-length_start     <= p.translation.y());
            REQUIRE(p.translation.y() <= length_start     );
            REQUIRE(-length_start     <= p.translation.z());
            REQUIRE(p.translation.z() <= length_start     );
            REQUIRE(-rad_start        <= p.rotation.x()   );
            REQUIRE(p.rotation.x()    <= rad_start        );
            REQUIRE(-rad_start        <= p.rotation.y()   );
            REQUIRE(p.rotation.y()    <= rad_start        );
            REQUIRE(-rad_start        <= p.rotation.z()   );
            REQUIRE(p.rotation.z()    <= rad_start        );
        }
    }
}