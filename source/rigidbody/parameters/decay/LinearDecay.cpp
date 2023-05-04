#include <rigidbody/parameters/decay/LinearDecay.h>

using namespace rigidbody::parameters::decay;

LinearDecay::LinearDecay(unsigned int max_iterations) {
    set_characteristic_time(max_iterations/2);
}

LinearDecay::LinearDecay(double decay_rate) : decay_rate(decay_rate) {}

double LinearDecay::get_factor() {
    return 1.0 - decay_rate*draws++;
}

void LinearDecay::set_characteristic_time(unsigned int iterations) {
    decay_rate = iterations/2.0;
}