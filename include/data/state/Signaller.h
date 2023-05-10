#pragma once

namespace signaller {
    /**
     * @brief A small probe for signalling changes which can be dispatched to other classes. 
     */
    class Signaller {
        public: 
            Signaller() {}

            /**
             * @brief Signal that the external state (i.e. position, rotation) of this object has changed. 
             */
            virtual void external_change() const = 0;

            /**
             * @brief Signal that the internal state (removed or added atoms) of this object has changed.
             */
            virtual void internal_change() const = 0;
    };
}
