#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/parameters/decay/ExponentialDecay.h>
#include <rigidbody/parameters/decay/LinearDecay.h>
#include <rigidbody/parameters/decay/NoDecay.h>
#include <rigidbody/Rigidbody.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("DecayFactory::create_decay_strategy") {
    SECTION("Linear") {
        auto decay = factory::create_decay_strategy(100, settings::rigidbody::DecayStrategyChoice::Linear);
        REQUIRE(dynamic_cast<parameter::decay::LinearDecay*>(decay.get()) != nullptr);
    }

    SECTION("Exponential") {
        auto decay = factory::create_decay_strategy(100, settings::rigidbody::DecayStrategyChoice::Exponential);
        REQUIRE(dynamic_cast<parameter::decay::ExponentialDecay*>(decay.get()) != nullptr);
    }

    SECTION("None") {
        auto decay = factory::create_decay_strategy(100, settings::rigidbody::DecayStrategyChoice::None);
        REQUIRE(dynamic_cast<parameter::decay::NoDecay*>(decay.get()) != nullptr);
    }
}
