// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/detail/InlineSignature.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/constraints/Constraint.h>
#include <utility/observer_ptr.h>

#include <vector>
#include <string>
#include <memory>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Set a relative level of hydration molecules for each body.
     */
    class RelativeHydrationElement : public GenericElement {
        public:
            /**
             * @brief Declare a hydration level for one body.
             */
            RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::string& name, double ratio);
            ~RelativeHydrationElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static InlineSignature _valid_inline_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

            /**
             * @brief The hydration weight of every body, indexed by body index, as handed to BodyCounterCulling.
             */
            std::vector<double> _get_ratios() const;

        private:
            observer_ptr<Sequencer> owner;
    };
}