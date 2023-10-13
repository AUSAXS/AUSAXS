#pragma once

#include <data/DataFwd.h>
#include <data/state/DataStateFwd.h>

#include <memory>

namespace signaller {class Signaller;}
namespace hist {
    struct BodyTracker {
        BodyTracker(const data::Molecule* const protein);

        ~BodyTracker();

        /**
         * @brief Get a signalling object for signalling a change of state. 
         *        Each body is supposed to hold one of these, and trigger it when they change state. 
         */
        std::shared_ptr<signaller::Signaller> get_probe(unsigned int i);

        /**
         * @brief Signal that the hydration layer was modified. 
         *        This is supposed to be used only by the Protein class, which has direct access to this object. Thus a signalling object is unnecessary. 
         */
        void signal_modified_hydration_layer();

        const state::StateManager* get_state_manager() const;

        state::StateManager* get_state_manager();

        const unsigned int body_size;                       // number of managed bodies
        std::unique_ptr<state::StateManager> statemanager;  // a helper which keeps track of state changes in each body
    };
}