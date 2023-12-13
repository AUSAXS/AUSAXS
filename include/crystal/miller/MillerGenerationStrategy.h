#pragma once

#include <crystal/miller/CrystalMillerFwd.h>

#include <vector>

namespace crystal {
    struct MillerGenerationStrategy {
        virtual std::vector<crystal::Miller> generate() const = 0;
    };
}