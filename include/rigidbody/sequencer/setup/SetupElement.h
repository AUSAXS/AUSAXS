#pragma once

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>

#include <string>
#include <vector>
#include <unordered_map>

namespace rigidbody::sequencer {
    /**
     * @brief Set up the optimization problem.
     *        Once any method from the LoopElementCallback is called, the setup is considered complete.
     */
    class SetupElement : public LoopElementCallback {
        public:
            SetupElement(observer_ptr<Sequencer> owner);
            virtual ~SetupElement() = default;

            SetupElement& load(const std::vector<std::string>& path, const std::vector<std::string>& body_names = {});

            SetupElement& load_existing(observer_ptr<RigidBody> rigidbody);

            SetupElement& distance_constraint();

            SetupElement& fixed_constraint();

            SetupElement& generate_linear_constraints();

            SetupElement& generate_volumetric_constraints();

            /**
             * @brief Get the name identifiers of all loaded bodies.
             */
            std::unordered_map<std::string, unsigned int>& _get_body_names();

            /**
             * @brief Set the currently active body for the setup.
             *        The currently active setup body will not influence the optimization in any way.
             */
            void _set_active_body(observer_ptr<RigidBody> body);

        private:
            std::unordered_map<std::string, unsigned int> body_names;
            observer_ptr<RigidBody> active_body;
    };
}