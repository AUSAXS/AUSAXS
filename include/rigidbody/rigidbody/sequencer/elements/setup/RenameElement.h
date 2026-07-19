// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    class RenameElement : public GenericElement {
        public:
            /**
             * @brief Rename an existing body.
             */
            RenameElement(observer_ptr<Sequencer> owner, std::string_view old_name, std::string_view new_name);

            ~RenameElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
    };
}
