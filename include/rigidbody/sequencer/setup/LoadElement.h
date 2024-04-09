#pragma once

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

#include <string>
#include <vector>
#include <memory>

namespace rigidbody::sequencer {
    class LoadElement : public GenericElement {
        public:
            LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names = {});
            LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& split, const std::vector<std::string>& body_names = {});
            ~LoadElement() override = default;

            void run() override;
        
        private:
            observer_ptr<Sequencer> owner;
            std::unique_ptr<rigidbody::RigidBody> rigidbody;
    };
}