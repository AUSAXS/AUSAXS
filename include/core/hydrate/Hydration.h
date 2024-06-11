#pragma once

namespace hydrate {
    class Hydration {
        public:
            Hydration() = default;
            virtual ~Hydration() = default;

            /**
             * @brief Clear the current hydration model.
             */
            virtual void clear() = 0;
    };
}