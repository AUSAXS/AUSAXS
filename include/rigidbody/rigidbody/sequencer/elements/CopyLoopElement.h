// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/elements/PhaseGuard.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief A duplicate of a previous loop.
     */
    class CopyLoopElement : public GenericElement, public RunPhaseGuard {
        public:
            CopyLoopElement(observer_ptr<LoopElement> owner, observer_ptr<LoopElement> target);
            ~CopyLoopElement() override;

            void run() override;

            observer_ptr<LoopElement> _get_target() const;

        private:
            observer_ptr<LoopElement> owner;
            observer_ptr<LoopElement> target;
    };
}