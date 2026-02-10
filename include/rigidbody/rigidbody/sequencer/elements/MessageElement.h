// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/elements/LoopElementCallback.h>

namespace ausaxs::rigidbody::sequencer {
    class MessageElement : public LoopElementCallback, public GenericElement {
        public:
            MessageElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, std::string_view message, bool log);
            ~MessageElement() override;

            void run() override;

        private:
            std::function<void()> message_func;
            std::function<std::string()> parse_user_msg(std::string_view msg) const;
    };
}