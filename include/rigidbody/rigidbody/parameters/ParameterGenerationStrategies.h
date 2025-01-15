#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace ausaxs::rigidbody::parameter {
    template<bool TRANSLATE, bool ROTATE, bool SYMMETRY>
    class LimitedParameterGenerator : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~LimitedParameterGenerator() override = default;

            Parameter next() override;
    };
}

template<bool TRANSLATE, bool ROTATE, bool SYMMETRY>
ausaxs::rigidbody::parameter::Parameter ausaxs::rigidbody::parameter::LimitedParameterGenerator<TRANSLATE, ROTATE, SYMMETRY>::next() {
    double scaling = decay_strategy->next();

    double dx = 0, dy = 0, dz = 0;
    if constexpr (TRANSLATE) {
        dx = translation_dist(generator)*scaling;
        dy = translation_dist(generator)*scaling;
        dz = translation_dist(generator)*scaling;
    }

    double dr1 = 0, dr2 = 0, dr3 = 0;
    if constexpr (ROTATE) {
        dr1 = rotation_dist(generator)*scaling;
        dr2 = rotation_dist(generator)*scaling;
        dr3 = rotation_dist(generator)*scaling;
    }

    return Parameter(Vector3(dx, dy, dz), dr1, dr2, dr3);
}

namespace ausaxs::rigidbody::parameter {
    using TranslationsOnly = LimitedParameterGenerator<true, false, false>;
    using RotationsOnly = LimitedParameterGenerator<false, true, false>;
    using SymmetryOnly = LimitedParameterGenerator<false, false, true>;
}