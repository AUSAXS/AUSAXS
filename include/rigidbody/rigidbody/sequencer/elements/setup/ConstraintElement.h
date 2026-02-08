// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/constraints/Constraint.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::rigidbody::sequencer {
    class ConstraintElement : public GenericElement {
        public:
            ConstraintElement(observer_ptr<Sequencer> owner, std::unique_ptr<rigidbody::constraints::Constraint> constraint);

            ~ConstraintElement() override = default;

            void run() override;
    };
}