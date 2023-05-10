#pragma once

#include <vector>

namespace crystal {
    class Miller;
    struct MillerGenerationStrategy {
        virtual std::vector<crystal::Miller> generate() const = 0;
    };
}