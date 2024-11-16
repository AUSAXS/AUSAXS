#pragma once

#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>

namespace ausaxs::rigidbody::sequencer {
    class OptimizeStepElement : public LoopElementCallback, public GenericElement {
        public:
            OptimizeStepElement(LoopElement* owner);
            ~OptimizeStepElement() override;

            void run() override;

            OptimizeStepElement& save_on_improvement(const io::File& path);

        private:
            std::vector<std::unique_ptr<GenericElement>> elements;
    };
}