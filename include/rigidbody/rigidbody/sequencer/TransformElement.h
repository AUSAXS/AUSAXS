#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    class TransformElement : public LoopElementCallback, public GenericElement {
        public:
            TransformElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::transform::TransformStrategy> strategy);
            ~TransformElement();

            void run() override;

        private:
            std::shared_ptr<rigidbody::transform::TransformStrategy> strategy;
    };
}