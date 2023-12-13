#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/parameters/decay/ExponentialDecay.h>
#include <rigidbody/parameters/decay/LinearDecay.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

std::unique_ptr<rigidbody::parameters::decay::DecayStrategy> rigidbody::factory::create_decay_strategy(unsigned int iterations) {
    return create_decay_strategy(iterations, settings::rigidbody::decay_strategy);
}

std::unique_ptr<rigidbody::parameters::decay::DecayStrategy> rigidbody::factory::create_decay_strategy(unsigned int iterations, const settings::rigidbody::DecayStrategyChoice& choice) {
    switch (choice) {
        case settings::rigidbody::DecayStrategyChoice::Linear:
            return std::make_unique<rigidbody::parameters::decay::LinearDecay>(iterations);
        case settings::rigidbody::DecayStrategyChoice::Exponential:
            return std::make_unique<rigidbody::parameters::decay::ExponentialDecay>(iterations);
        default:
            throw except::unknown_argument("rigidbody::parameters::factory::create_decay_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}