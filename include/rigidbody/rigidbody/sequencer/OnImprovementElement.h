#pragma once

#include <rigidbody/sequencer/LoopElement.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    class OnImprovementElement : public LoopElement {
        public:
            OnImprovementElement(observer_ptr<LoopElement> owner);
            ~OnImprovementElement() override;

            void run() override;

        private:
            double best_chi2;
    };
}