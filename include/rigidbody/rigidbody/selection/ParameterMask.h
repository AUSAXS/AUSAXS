// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/BodyTransformParametersRelative.h>

namespace ausaxs::rigidbody::selection {
    /**
     * @brief A mask that filters generated transformation parameters, keeping only selected components active.
     *
     * Call apply() after generating parameters and before applying them to the body.
     * The apply() method asserts that at least one parameter field remains active after masking,
     * catching mismatches between the mask configuration and the parameter generator.
     */
    struct ParameterMask {
        bool real_translation  = true; // Allow body translation and rotation.
        bool real_rotation     = true; // Allow body rotation.
        bool sym_translation   = true; // Allow symmetry offset translation (initial_relation.translation).
        bool sym_axis          = true; // Allow symmetry rotation axis direction (repeat_relation.axis).

        static ParameterMask all()                  { return {true,  true,  true,  true }; }
        static ParameterMask real_only()            { return {true,  true,  false, false}; }
        static ParameterMask symmetry_only()        { return {false, false, true,  true }; }
        static ParameterMask real_only_rot ()       { return {false, true,  false, false}; }
        static ParameterMask real_only_trans ()     { return {true,  false, false, false}; }
        static ParameterMask symmetry_only_trans()  { return {false, false, true,  false}; }
        static ParameterMask symmetry_only_axis()   { return {false, false, false, true }; }

        /**
         * @brief Apply this mask to relative transform parameters in-place.
         *
         * Masked-out fields are cleared: optional fields become std::nullopt,
         * symmetry sub-fields are set to zero delta. Asserts that at least one
         * parameter field remains active after masking.
         */
        void apply(parameter::BodyTransformParametersRelative& params) const;
    };
}
