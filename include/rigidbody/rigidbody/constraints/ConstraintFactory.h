// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/Constraint.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<constraints::Constraint> create_constraint_cm(
        observer_ptr<const data::Molecule> owner, sequencer::detail::BodySymmetrySelector body1, sequencer::detail::BodySymmetrySelector body2
    );

    std::unique_ptr<constraints::Constraint> create_constraint_bond(
        observer_ptr<const data::Molecule> owner, sequencer::detail::BodySymmetrySelector body1, sequencer::detail::BodySymmetrySelector body2
    );

    std::unique_ptr<constraints::Constraint> create_constraint_attractor(
        observer_ptr<const data::Molecule> owner, sequencer::detail::BodySymmetrySelector body1, sequencer::detail::BodySymmetrySelector body2, double target_distance
    );

    std::unique_ptr<constraints::Constraint> create_constraint_repeller(
        observer_ptr<const data::Molecule> owner, sequencer::detail::BodySymmetrySelector body1, sequencer::detail::BodySymmetrySelector body2, double target_distance
    );

    std::unique_ptr<constraints::Constraint> create_constraint(
        observer_ptr<const data::Molecule> owner, sequencer::detail::BodySymmetrySelector body1, sequencer::detail::BodySymmetrySelector body2, unsigned int iatom1, unsigned int iatom2
    );
}