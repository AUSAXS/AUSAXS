// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <mini/Minimizer.h>

#include <vector>

namespace ausaxs::mini {
    /**
     * @brief A Levenberg-Marquardt least-squares minimizer.
     *
     * The Jacobian of the residuals is approximated by forward differences, so each iteration costs one evaluation per free parameter
     * plus one per trial step. Bounded parameters are handled by projecting each step onto the box, and freezing any parameter which
     * sits on a bound while the gradient points out of it.
     */
    class LevenbergMarquardt : public Minimizer {
        public:
            LevenbergMarquardt() = default;
            LevenbergMarquardt(residual_function func, const std::vector<Parameter>& params);
            ~LevenbergMarquardt() override = default;

        private:
            Result minimize_override() override;
    };
}
