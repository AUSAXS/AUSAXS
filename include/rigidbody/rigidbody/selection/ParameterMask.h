// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/BodyTransformParametersRelative.h>

#include <optional>

namespace ausaxs::rigidbody::selection {
    /**
     * @brief A mask that filters generated transformation parameters, keeping only selected components active.
     *
     * Call apply() after generating parameters and before applying them to the body to assert at least one parameter field is active after masking. 
     */
    struct ParameterMask {
        bool real_translation  = true; // Allow body translation and rotation.
        bool real_rotation     = true; // Allow body rotation.
        bool sym_translation   = true; // Allow symmetry offset translation (initial_relation.translation).
        bool sym_axis          = true; // Allow symmetry rotation axis direction (repeat_relation.axis).

        // If set, the two symmetry flags above only apply to the body's symmetry at this index. Unset means they apply to all symmetries alike.
        std::optional<int> target_symmetry = std::nullopt;

        static ParameterMask all()                  { return {.real_translation=true,  .real_rotation=true,  .sym_translation=true,  .sym_axis=true }; }
        static ParameterMask real_only()            { return {.real_translation=true,  .real_rotation=true,  .sym_translation=false, .sym_axis=false}; }
        static ParameterMask symmetry_only()        { return {.real_translation=false, .real_rotation=false, .sym_translation=true,  .sym_axis=true }; }
        static ParameterMask real_only_rot ()       { return {.real_translation=false, .real_rotation=true,  .sym_translation=false, .sym_axis=false}; }
        static ParameterMask real_only_trans ()     { return {.real_translation=true,  .real_rotation=false, .sym_translation=false, .sym_axis=false}; }
        static ParameterMask symmetry_only_trans()  { return {.real_translation=false, .real_rotation=false, .sym_translation=true,  .sym_axis=false}; }
        static ParameterMask symmetry_only_axis()   { return {.real_translation=false, .real_rotation=false, .sym_translation=false, .sym_axis=true }; }

        /**
         * @brief Apply this mask to relative transform parameters in-place.
         *
         * Masked-out fields are cleared: optional fields become std::nullopt, symmetry sub-fields are set to zero delta. 
         */
        void apply(parameter::BodyTransformParametersRelative& params) const;
    };
}
