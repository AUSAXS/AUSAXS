#pragma once

class Protein;
class StateManager;
namespace hist {
    struct BodyTracker {
        BodyTracker(Protein* protein);

        ~BodyTracker();

        /**
         * @brief Get a signalling object for signalling a change of state. 
         *        Each body is supposed to hold one of these, and trigger it when they change state. 
         */
        std::shared_ptr<StateManager::BoundSignaller> get_probe(unsigned int i);

        /**
         * @brief Signal that the hydration layer was modified. 
         *        This is supposed to be used only by the Protein class, which has direct access to this object. Thus a signalling object is unnecessary. 
         */
        void signal_modified_hydration_layer();

        const StateManager& get_state_manager() const;

        StateManager& get_state_manager();

        const unsigned int size;                    // number of managed bodies
        std::unique_ptr<StateManager> statemanager; // a helper which keeps track of state changes in each body
    };
}