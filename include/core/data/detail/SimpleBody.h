#pragma once

#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>

namespace ausaxs::data::detail {
    struct SimpleBody {
        std::vector<AtomFF> atoms;
        std::vector<Water> waters;
    };
}