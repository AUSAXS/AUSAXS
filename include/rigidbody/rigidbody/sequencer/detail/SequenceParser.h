// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/Sequencer.h>

namespace ausaxs::rigidbody::sequencer {
    class SequenceParser {
        public:
            SequenceParser() = default;

            std::unique_ptr<Sequencer> parse_file(const io::ExistingFile& config);
            std::unique_ptr<Sequencer> parse_text(const std::string& script);

        private:
            std::unique_ptr<Sequencer> parse(std::istream& in, const std::string& config_folder = "");
            std::vector<observer_ptr<LoopElement>> loop_stack;
    };
}