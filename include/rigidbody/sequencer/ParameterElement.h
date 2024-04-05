#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>

namespace rigidbody::sequencer {
    class ParameterElement : public LoopElementCallback, public GenericElement {
        public:
            ParameterElement(LoopElement* owner, std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy);
            ~ParameterElement() override;

            void run() override;

            ParameterElement& max_rotation_angle(double radians);

            ParameterElement& max_translation_distance(double distance);

            ParameterElement& decay_strategy(std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> strategy);

        private: 
            std::shared_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy;
    };
}