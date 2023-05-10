#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/parameters/decay/ExponentialDecay.h>
#include <rigidbody/parameters/decay/LinearDecay.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace rigidbody::parameters::factory;

std::unique_ptr<rigidbody::parameters::decay::DecayStrategy> rigidbody::parameters::factory::create_decay_strategy() {
    return create_decay_strategy(settings::rigidbody::decay_strategy);
}

std::unique_ptr<rigidbody::parameters::decay::DecayStrategy> rigidbody::parameters::factory::create_decay_strategy(const settings::rigidbody::DecayStrategyChoice& choice) {
    switch (choice) {
        case settings::rigidbody::DecayStrategyChoice::Linear:
            return std::make_unique<rigidbody::parameters::decay::LinearDecay>();
        case settings::rigidbody::DecayStrategyChoice::Exponential:
            return std::make_unique<rigidbody::parameters::decay::ExponentialDecay>();
        default:
            throw except::unknown_argument("rigidbody::parameters::factory::create_decay_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}