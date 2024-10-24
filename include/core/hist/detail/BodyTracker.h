#pragma once

#include <data/DataFwd.h>
#include <data/state/DataStateFwd.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::signaller {class Signaller;}
namespace ausaxs::hist {
    struct BodyTracker {
        BodyTracker(observer_ptr<const data::Molecule> protein);

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

        observer_ptr<const state::StateManager> get_state_manager() const;

        observer_ptr<state::StateManager> get_state_manager();

        const unsigned int body_size;                       // number of managed bodies
        std::unique_ptr<state::StateManager> statemanager;  // a helper which keeps track of state changes in each body
    };
}