// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <crystal/miller/CrystalMillerFwd.h>

#include <vector>

namespace ausaxs::crystal {
    struct MillerGenerationStrategy {
        virtual ~MillerGenerationStrategy() = default;
        virtual std::vector<crystal::Miller> generate() const = 0;
    };
}