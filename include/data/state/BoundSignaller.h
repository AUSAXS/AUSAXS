#pragma once

#include <data/state/Signaller.h>

class StateManager;
namespace signaller {
    /**
     * @brief A small probe for signalling changes which can be dispatched to other classes. 
     */
    class BoundSignaller : public Signaller {
        public: 
            BoundSignaller(unsigned int id, StateManager* const owner);

            ~BoundSignaller() override;

            /**
             * @brief Signal that the external state (i.e. position, rotation) of this object has changed. 
             */
            virtual void external_change() const;

            /**
             * @brief Signal that the internal state (removed or added atoms) of this object has changed.
             */
            virtual void internal_change() const;

            /**
             * @brief Get the id of this signaller. 
             */
            unsigned int get_id() const;

        private: 
            StateManager* const owner;
            unsigned int id;
    };
}