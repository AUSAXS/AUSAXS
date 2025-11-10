// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/xray/FormFactor.h>
#include <form_factor/xray/ExvFormFactor.h>
#include <form_factor/tables/FormFactorTableXray.h>
#include <constants/ConstantsFwd.h>
#include <data/DataFwd.h>

#include <array>
#include <stdexcept>

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
            constexpr FormFactor(ExvFormFactor&& ffx) : a({ffx.q0, 0, 0, 0, 0}), b({ffx.exponent, 0, 0, 0, 0}), c(0), f0(1) {}

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
                this->f0 = f0/I0();
            }

        private: 
            std::array<double, 5> a;
            std::array<double, 5> b;
            double c;
            double f0 = 1;
    };

    /**
     * This struct contains the form factors of the most common atomic elements encountered in SAXS. 
     */
    namespace storage::atomic {
        using namespace ausaxs::constants::form_factor;
        // atomic
        constexpr FormFactor H               = FormFactor(               xray::H::a,               xray::H::b,               xray::H::c);
        constexpr FormFactor C               = FormFactor(               xray::C::a,               xray::C::b,               xray::C::c);
        constexpr FormFactor N               = FormFactor(               xray::N::a,               xray::N::b,               xray::N::c);
        constexpr FormFactor O               = FormFactor(               xray::O::a,               xray::O::b,               xray::O::c);
        constexpr FormFactor S               = FormFactor(               xray::S::a,               xray::S::b,               xray::S::c);

        // atomic groups
        constexpr FormFactor CH_sp3          = FormFactor(          xray::CH_sp3::a,          xray::CH_sp3::b,          xray::CH_sp3::c);
        constexpr FormFactor CH2_sp3         = FormFactor(         xray::CH2_sp3::a,         xray::CH2_sp3::b,         xray::CH2_sp3::c);
        constexpr FormFactor CH3_sp3         = FormFactor(         xray::CH3_sp3::a,         xray::CH3_sp3::b,         xray::CH3_sp3::c);
        constexpr FormFactor CH_sp2          = FormFactor(          xray::CH_sp2::a,          xray::CH_sp2::b,          xray::CH_sp2::c);
        constexpr FormFactor CH_arom         = FormFactor(         xray::CH_arom::a,         xray::CH_arom::b,         xray::CH_arom::c);
        constexpr FormFactor OH_alc          = FormFactor(          xray::OH_alc::a,          xray::OH_alc::b,          xray::OH_alc::c);
        constexpr FormFactor OH_acid         = FormFactor(         xray::OH_acid::a,         xray::OH_acid::b,         xray::OH_acid::c);
        constexpr FormFactor O_res           = FormFactor(           xray::O_res::a,           xray::O_res::b,           xray::O_res::c);
        constexpr FormFactor NH              = FormFactor(              xray::NH::a,              xray::NH::b,              xray::NH::c);
        constexpr FormFactor NH2             = FormFactor(             xray::NH2::a,             xray::NH2::b,             xray::NH2::c);
        constexpr FormFactor NH_plus         = FormFactor(         xray::NH_plus::a,         xray::NH_plus::b,         xray::NH_plus::c);
        constexpr FormFactor NH2_plus        = FormFactor(        xray::NH2_plus::a,        xray::NH2_plus::b,        xray::NH2_plus::c);
        constexpr FormFactor NH3_plus        = FormFactor(        xray::NH3_plus::a,        xray::NH3_plus::b,        xray::NH3_plus::c);
        constexpr FormFactor NH_guanine      = FormFactor(      xray::NH_guanine::a,      xray::NH_guanine::b,      xray::NH_guanine::c);
        constexpr FormFactor NH2_guanine     = FormFactor(     xray::NH2_guanine::a,     xray::NH2_guanine::b,     xray::NH2_guanine::c);
        constexpr FormFactor SH              = FormFactor(              xray::SH::a,              xray::SH::b,              xray::SH::c);

        // average excluded volume
        constexpr FormFactor excluded_volume = FormFactor( xray::excluded_volume::a, xray::excluded_volume::b, xray::excluded_volume::c);

        // all others; this is just the form factor of argon
        constexpr FormFactor other           = FormFactor(           xray::other::a,           xray::other::b,           xray::other::c);

        constexpr const FormFactor& get_form_factor(form_factor_t type) {
            switch (type) {
                case form_factor_t::H:                  return H;
                case form_factor_t::C:                  return C;
                case form_factor_t::N:                  return N;
                case form_factor_t::O:                  return O;
                case form_factor_t::S:                  return S;
                case form_factor_t::CH:                 return CH_sp3;
                case form_factor_t::CH2:                return CH2_sp3;
                case form_factor_t::CH3:                return CH3_sp3;
                case form_factor_t::NH:                 return NH;
                case form_factor_t::NH2:                return NH2;
                case form_factor_t::NH3:                return NH3_plus;
                case form_factor_t::OH:                 return OH_alc;
                case form_factor_t::SH:                 return SH;
                case form_factor_t::OTHER:              return other;
                case form_factor_t::EXCLUDED_VOLUME:    return excluded_volume;
                default: throw std::runtime_error("form_factor::storage::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
        }
    }
}