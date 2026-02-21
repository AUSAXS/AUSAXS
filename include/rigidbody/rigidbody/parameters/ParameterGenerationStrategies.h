// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace ausaxs::rigidbody::parameter {
    template<bool TRANSLATE, bool ROTATE, bool SYMMETRY>
    class LimitedParameterGenerator : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~LimitedParameterGenerator() override = default;

            BodyTransformParametersRelative next(int ibody) override;
    };
}

namespace ausaxs::rigidbody::parameter {
    using AllParameters    = LimitedParameterGenerator<true,  true,  true>;
    using TranslationsOnly = LimitedParameterGenerator<true,  false, false>;
    using RotationsOnly    = LimitedParameterGenerator<false, true,  false>;
    using SymmetryOnly     = LimitedParameterGenerator<false, false, true>;
}