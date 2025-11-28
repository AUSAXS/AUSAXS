// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/ExvFormFactor.h>
#include <constants/ConstantsFwd.h>
#include <data/DataFwd.h>

#include <array>

namespace ausaxs::form_factor {
    class FormFactor {
        public:
            /**
             * @brief Initialize a vacuum form factor based on a 5-Gaussian approximation.
             */
            constexpr FormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {}

            /**
             * @brief Initialize an excluded volume form factor.
             *        This is only used to instantiate the average excluded volume form factor.
             *        Note that these excluded volume form factors are not normalized. 
             */
            constexpr FormFactor(ExvFormFactor&& ffx) : a({ffx.q0, 0, 0, 0, 0}), b({ffx.exponent, 0, 0, 0, 0}), c(0) {}

            /**
             * @brief Evaluate the form factor at a given q value.
             *        The vacuum form factors are normalized to 1 at q = 0.
             */
            constexpr double evaluate(double q) const {
                double sum = 0;
                for (unsigned int i = 0; i < 5; ++i) {
                    sum += a[i]*constexpr_math::exp(-b[i]*q*q);
                }
                return (sum + c)*f0;
            }

            /**
             * @brief Evaluate the form factor at q = 0.
             */
            constexpr double I0() const {
                return (a[0] + a[1] + a[2] + a[3] + a[4] + c);
            }

            /**
             * @brief Manually set the normalization of this form factor.
             *        evaluate(0) will return this value.
             */
            constexpr void set_normalization(double f0) {
                this->f0 = f0;
            }

        protected:
            double f0 = 1;

        private: 
            std::array<double, 5> a;
            std::array<double, 5> b;
            double c;
    };
}