// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    class BodySelectElement : public LoopElementCallback, public GenericElement {
        public:
            BodySelectElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy);
            ~BodySelectElement();

            void run() override;

        private:
            std::shared_ptr<rigidbody::selection::BodySelectStrategy> strategy;
    };
}