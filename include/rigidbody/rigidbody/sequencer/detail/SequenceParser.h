// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/Sequencer.h>

#include <unordered_map>

namespace ausaxs::rigidbody::sequencer {
    enum class ElementType;

    class SequenceParser {
        public:
            SequenceParser() = default;

            std::unique_ptr<Sequencer> parse(const io::ExistingFile& config);

        private:
            std::vector<observer_ptr<LoopElement>> loop_stack;

            template<ElementType T>
            std::unique_ptr<GenericElement> parse_arguments(const std::unordered_map<std::string, std::vector<std::string>>& args);
    };
}