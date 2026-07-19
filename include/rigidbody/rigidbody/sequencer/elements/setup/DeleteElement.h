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
             *        This is a setup-time operation: it takes effect immediately when parsed/constructed, so it
             *        must appear before any element that refers to bodies by index or name (e.g. symmetry,
             *        constraints). Doing otherwise is a user error and will lead to invalid references.
             */
            DeleteElement(observer_ptr<Sequencer> owner, std::vector<std::string> names);

            ~DeleteElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);
    };
}
