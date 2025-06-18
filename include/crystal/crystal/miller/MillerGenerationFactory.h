// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>
#include <settings/CrystalSettings.h>

#include <memory>

namespace ausaxs::crystal {
    namespace factory {
        /**
         * @brief Prepare a miller generation strategy.
         */
        std::unique_ptr<crystal::MillerGenerationStrategy> construct_miller_strategy();

        /**
         * @brief Prepare a miller generation strategy.
         */
        std::unique_ptr<crystal::MillerGenerationStrategy> construct_miller_strategy(const settings::crystal::MillerGenerationChoice& choice);
    }
}