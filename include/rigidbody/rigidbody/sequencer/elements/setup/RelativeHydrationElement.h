// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/constraints/Constraint.h>
#include <utility/observer_ptr.h>

#include <vector>
#include <string>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Set a relative level of hydration molecules for each body.
     */
    class RelativeHydrationElement : public GenericElement {
        public:
            RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<double>& ratios);
            ~RelativeHydrationElement() override;

            void run() override;

        private:
            observer_ptr<Sequencer> owner;
            std::vector<double> ratios;
    };
}