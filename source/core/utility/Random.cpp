#include <utility/Random.h>

using namespace ausaxs;

std::mt19937 gen{std::random_device{}()};

void random::set_seed(int seed) noexcept {
    gen.seed(seed);
}

const std::mt19937& random::generator() {
    return gen;
}