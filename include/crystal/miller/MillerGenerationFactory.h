#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>

#include <memory>

namespace settings::crystal {enum class MillerGenerationChoice;}
namespace crystal {
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