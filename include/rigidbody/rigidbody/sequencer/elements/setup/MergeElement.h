// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/detail/InlineSignature.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    class MergeElement : public GenericElement {
        public:
            /**
             * @brief Merge the atoms of one or more existing bodies into another, deleting the merged-away bodies.
             */
            MergeElement(observer_ptr<Sequencer> owner, std::string_view first_name, std::vector<std::string> other_names);

            ~MergeElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static InlineSignature _valid_inline_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
    };
}
