#pragma once

#include <data/state/DataStateFwd.h>
#include <data/state/Signaller.h>

namespace signaller {
    /**
     * @brief A small probe for signalling changes which can be dispatched to other classes. 
     */
    class BoundSignaller : public Signaller {
        public: 
            BoundSignaller(unsigned int id, state::StateManager* const owner);

            ~BoundSignaller() override;

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