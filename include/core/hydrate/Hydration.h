#pragma once

#include <data/atoms/Water.h>

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

            /**
             * @brief Create a new hydration data structure. The type will depend on the argument.
             */
            template<data::WaterVector T>
            static std::unique_ptr<Hydration> create(T&& data);

            static std::unique_ptr<Hydration> create(); //< @copydoc create()
    };
}