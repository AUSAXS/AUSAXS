#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>
#include <settings/CrystalSettings.h>

namespace crystal {
    namespace factory {
        /**
         * @brief Prepare a miller generation strategy.
         */
        std::unique_ptr<crystal::MillerGenerationStrategy> construct_miller_strategy(settings::crystal::MillerGenerationChoice choice = settings::crystal::miller_generation_strategy);
    }
}