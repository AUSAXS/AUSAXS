// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    enum class ElementType {
        AutomaticConstraint,
        BodySelect,
        Constraint,
        Copy,
        EveryNStep,
        LoadElement,
        Log,
        LoopBegin,
        LoopEnd,
        Message,
        OnImprovement,
        OptimizeStep,
        OutputFolder,
        OverlapStrength,
        Parameter,
        RelativeHydration,
        Save,
        Seed,
        SymmetryElement,
        Transform,
        COUNT
    };

    std::vector<std::string> valid_elements();
    std::vector<std::string> valid_arguments(ElementType type);
    ElementType get_type(std::string_view line);
}