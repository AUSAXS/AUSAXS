#pragma once

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>

namespace settings::rigidbody {enum class ConstraintGenerationStrategyChoice;}

namespace rigidbody::sequencer {
    class AutoConstraintsElement : public GenericElement {
        public:
            AutoConstraintsElement(observer_ptr<Sequencer> owner, settings::rigidbody::ConstraintGenerationStrategyChoice strategy);
            ~AutoConstraintsElement() override = default;

            void run() override;

        private: 
            observer_ptr<Sequencer> owner;
            settings::rigidbody::ConstraintGenerationStrategyChoice strategy;
    };
}