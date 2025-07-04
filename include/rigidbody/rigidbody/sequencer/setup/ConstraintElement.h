// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/constraints/Constraint.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>

namespace ausaxs::rigidbody::sequencer {
    class ConstraintElement : public GenericElement {
        public:
            ConstraintElement(observer_ptr<Sequencer> owner, std::unique_ptr<rigidbody::constraints::Constraint> constraint);
            ConstraintElement(observer_ptr<Sequencer> owner, const std::string& body1, const std::string& body2, bool center_mass = false);
            ConstraintElement(observer_ptr<Sequencer> owner, const std::string& body1, const std::string& body2, unsigned int iatom1, unsigned int iatom2);

            ~ConstraintElement() override = default;

            void run() override;
    };
}