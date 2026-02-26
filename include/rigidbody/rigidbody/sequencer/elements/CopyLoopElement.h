// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief A duplicate of a previous loop.
     */
    class CopyLoopElement : public GenericElement {
        public:
            CopyLoopElement(observer_ptr<LoopElement> owner, observer_ptr<LoopElement> target);
            virtual ~CopyLoopElement();

            void run() override;

        private:
            observer_ptr<LoopElement> owner;
            observer_ptr<LoopElement> target;
    };
}