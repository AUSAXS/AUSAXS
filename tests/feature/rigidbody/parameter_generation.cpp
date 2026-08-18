#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

#include <algorithm>
#include <span>

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
TEST_CASE("UniformParameterGenerator::next symmetry components are independent") {
    settings::general::verbose = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    int iterations = 100;
    Molecule m{std::vector<Body>{Body(std::vector{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C), AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
    })}};
    m.get_body(0).symmetry().add(symmetry::type::c2);
    Rigidbody rb(std::move(m));

    auto all_zero = [](std::span<double> s) {return std::all_of(s.begin(), s.end(), [](double v) {return v == 0;});};
    auto any_nonzero = [](std::span<double> s) {return std::any_of(s.begin(), s.end(), [](double v) {return v != 0;});};

    SECTION("only the symmetry translation moves when only its amplitude is set") {
        rigidbody::parameter::UniformParameterGenerator gen(&rb, iterations, {.symmetry_translation = 5});

        bool saw_translation = false;
        for (int i = 0; i < iterations; i++) {
            auto p = gen.next(0);
            CHECK_FALSE(p.translation.has_value());
            CHECK_FALSE(p.rotation.has_value());
            REQUIRE(p.symmetry_pars.has_value());
            REQUIRE(p.symmetry_pars.value().size() == 1);
            auto& delta = *p.symmetry_pars.value()[0];
            CHECK(all_zero(delta.span_rotation()));
            saw_translation |= any_nonzero(delta.span_translation());
        }
        CHECK(saw_translation);
    }

    SECTION("only the symmetry rotation moves when only its amplitude is set") {
        rigidbody::parameter::UniformParameterGenerator gen(&rb, iterations, {.symmetry_rotation = 0.5});

        bool saw_rotation = false;
        for (int i = 0; i < iterations; i++) {
            auto p = gen.next(0);
            CHECK_FALSE(p.translation.has_value());
            CHECK_FALSE(p.rotation.has_value());
            REQUIRE(p.symmetry_pars.has_value());
            REQUIRE(p.symmetry_pars.value().size() == 1);
            auto& delta = *p.symmetry_pars.value()[0];
            CHECK(all_zero(delta.span_translation()));
            saw_rotation |= any_nonzero(delta.span_rotation());
        }
        CHECK(saw_rotation);
    }

    SECTION("no symmetry parameters are generated when neither amplitude is set") {
        rigidbody::parameter::UniformParameterGenerator gen(&rb, iterations, {.translation = 1, .rotation = 0.1});
        for (int i = 0; i < iterations; i++) {
            auto p = gen.next(0);
            CHECK(p.translation.has_value());
            CHECK(p.rotation.has_value());
            CHECK_FALSE(p.symmetry_pars.has_value());
        }
    }
}
