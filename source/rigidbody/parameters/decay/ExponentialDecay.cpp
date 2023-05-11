#include <rigidbody/parameters/decay/ExponentialDecay.h>

#include <cmath>

using namespace rigidbody::parameters::decay;

ExponentialDecay::ExponentialDecay(unsigned int max_iterations) {
    set_characteristic_time(max_iterations/2);
}

ExponentialDecay::ExponentialDecay(double decay_rate) : decay_rate(decay_rate) {}

ExponentialDecay::~ExponentialDecay() = default;

double ExponentialDecay::get_factor() {
    return std::exp(-decay_rate*draws++);
}

void ExponentialDecay::set_characteristic_time(unsigned int iterations) {
    decay_rate = 1.0/iterations;
}