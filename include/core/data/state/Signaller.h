#pragma once

#include <cstddef>

namespace ausaxs::signaller {
    /**
     * @brief A small probe for signalling changes which can be dispatched to other classes. 
     */
    class Signaller {
        public: 
            Signaller() = default;
            virtual ~Signaller() = default;

            /**
             * @brief Signal that the external state (i.e. position, rotation) of this object has changed. 
             */
            virtual void modified_external() const = 0;

            /**
             * @brief Signal that the internal state (removed or added atoms) of this object has changed.
             */
            virtual void modified_internal() const = 0;

            /**
             * @brief Signal that the ith symmetry of this object has changed. 
             */
            virtual void modified_symmetry(int i) const = 0;

            /**
             * @brief Set the number of symmetries to track. 
             */
            virtual void set_symmetry_size(std::size_t size) const = 0;
    };
}
