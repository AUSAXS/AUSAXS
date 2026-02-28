// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/Sequencer.h>

namespace ausaxs::rigidbody::sequencer {
    class SequenceParser {
        public:
            SequenceParser() = default;

            std::unique_ptr<Sequencer> parse(const io::ExistingFile& config);

        private:
            std::vector<observer_ptr<LoopElement>> loop_stack;
    };
}