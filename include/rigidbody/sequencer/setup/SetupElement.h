#pragma once

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>

#include <vector>
#include <string>

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

            SetupElement& constraint();

        private:
            std::vector<std::string> body_names;
    };
}