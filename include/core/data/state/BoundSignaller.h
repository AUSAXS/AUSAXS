#pragma once

#include <data/state/DataStateFwd.h>
#include <data/state/Signaller.h>

namespace ausaxs::signaller {
    /**
     * @brief A small probe for signalling changes which can be dispatched to other classes. 
     */
    class BoundSignaller : public Signaller {
        public: 
			BoundSignaller(const BoundSignaller& rhs) = default;
			BoundSignaller(BoundSignaller&& rhs) noexcept = default;
			~BoundSignaller() override = default;

            BoundSignaller(unsigned int id, state::StateManager* const owner);

            /**
             * @brief Signal that the external state (i.e. position, rotation) of this object has changed. 
             */
            virtual void external_change() const override;

            /**
             * @brief Signal that the internal state (removed or added atoms) of this object has changed.
             */
            virtual void internal_change() const override;

            /**
             * @brief Get the id of this signaller. 
             */
            unsigned int get_id() const;

        private: 
            state::StateManager* const owner;
            unsigned int id;
    };
}