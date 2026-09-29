// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/detail/InlineSignature.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <utility/observer_ptr.h>

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    class LoopElement;

    class MessageElement : public GenericElement {
        public:
            MessageElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, std::string_view message, std::string_view colour, bool log);
            MessageElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, std::string_view message, bool log);
            ~MessageElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static InlineSignature _valid_inline_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            observer_ptr<LoopElement> owner;
            std::function<void()> message_func;
            std::function<std::string()> parse_user_msg(std::string_view msg) const;
    };
}