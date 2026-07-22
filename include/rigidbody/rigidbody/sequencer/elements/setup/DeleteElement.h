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
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    class DeleteElement : public GenericElement {
        public:
            /**
             * @brief Delete one or more existing bodies.
             */
            DeleteElement(observer_ptr<Sequencer> owner, std::vector<std::string> names);

            ~DeleteElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
    };
}
