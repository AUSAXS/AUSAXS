// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>
#include <io/ExistingFile.h>

#include <memory>

namespace ausaxs::rigidbody::sequencer {
    class CopyBodyElement : public GenericElement {
        public:
            /**
             * @brief Create a copy of an existing body. The new body will be shifted to avoid overlapping. 
             */
            CopyBodyElement(observer_ptr<Sequencer> owner, std::string_view body_name, std::string_view source_body_name);

            /**
             * @brief Create a copy of an existing body.
             */
            CopyBodyElement(observer_ptr<Sequencer> owner, std::string_view body_name, int source_body_index);

            ~CopyBodyElement() override;

            void run() override;

            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        private:
            observer_ptr<Sequencer> owner;
    };
}