#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Convert all symmetries of the given body to real bodies at their current positions.
     *        This will allow them to be further optimized as independent bodies. 
     */
    class ConvertSymmetryElement : public LoopElementCallback, public GenericElement {
        public:
            ConvertSymmetryElement(observer_ptr<LoopElement> owner, int body);
            ~ConvertSymmetryElement() override;

            void run() override;
    };
}