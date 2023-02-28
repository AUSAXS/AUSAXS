#pragma once

#include <crystal/miller/Miller.h>

#include <vector>

namespace crystal {
    struct MillerGenerationStrategy {
        virtual std::vector<Miller<>> generate() = 0;
    };
}