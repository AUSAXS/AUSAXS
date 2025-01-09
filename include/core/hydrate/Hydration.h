#pragma once

#include <memory>

namespace ausaxs::hydrate {
    class Hydration {
        public:
            Hydration() = default;
            virtual ~Hydration() = default;

            /**
             * @brief Clear the current hydration model.
             */
            virtual void clear() = 0;

            /**
             * @brief Clone this strategy. 
             */
            virtual std::unique_ptr<Hydration> clone() const = 0;
    };
}