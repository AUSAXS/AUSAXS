// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/detail/InlineSignature.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    class SettingsElement : public GenericElement {
        public:
            SettingsElement(observer_ptr<Sequencer> owner, std::string name, std::string value);
            ~SettingsElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static InlineSignature _valid_inline_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            observer_ptr<Sequencer> owner;
            std::string name;
            std::string value;
    };
}
