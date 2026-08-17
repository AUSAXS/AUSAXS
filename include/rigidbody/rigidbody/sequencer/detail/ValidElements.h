// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/detail/ParsedArgs.h>

#include <string>
#include <string_view>
#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    enum class ElementType {
        AutomaticConstraint,
        BodySelect,
        Constraint,
        ConvertToSymmetry,
        Copy,
        Delete,
        EveryNStep,
        LoadElement,
        Log,
        LoopBegin,
        LoopEnd,
        Merge,
        Message,
        OnImprovement,
        OptimizeStep,
        OutputFolder,
        OverlapStrength,
        Parameter,
        RelativeHydration,
        Rename,
        Save,
        Seed,
        Split,
        SymmetryElement,
        Transform,
        Update,
        COUNT
    };

    std::vector<std::string> valid_elements();
    std::vector<std::string> valid_arguments(ElementType type);
    ElementType get_type(std::string_view line);

    /**
     * @brief Reject every named argument whose key the element does not accept.
     * @param type The element type the arguments were parsed for.
     * @param element The element name as written in the script; used for the error message.
     */
    void validate_named_arguments(ElementType type, std::string_view element, const ParsedArgs& args);
}