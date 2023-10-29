#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactor.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/DisplacedVolumeTable.h>
#include <form_factor/ExvFormFactor.h>
#include <constants/ConstantsFwd.h>
#include <data/DataFwd.h>

#include <array>
#include <string_view>
#include <stdexcept>
#include <cmath>

namespace form_factor {
    class FormFactor {
        public:
            /**
             * @brief Initialize a vacuum form factor based on a 5-Gaussian approximation.
             */
            constexpr FormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {
                f0 = 1./(a[0] + a[1] + a[2] + a[3] + a[4] + c);
            }

            /**
             * @brief Initialize an excluded volume form factor.
             *        This is only used to instantiate the average excluded volume form factor.
             */
            constexpr FormFactor(ExvFormFactor&& ffx) : a({1, 0, 0, 0, 0}), b({ffx.exponent, 0, 0, 0, 0}), c(0), f0(1) {}

            /**
             * @brief Evaluate the form factor at a given q value.
             *        The vacuum form factors are normalized to 1 at q = 0.
             */
            constexpr double evaluate(double q) const {
                double sum = 0;
                for (unsigned int i = 0; i < 5; ++i) {
                    sum += a[i]*std::exp(-b[i]*q*q);
                }
                return (sum + c)*f0;
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
        // atomic
        constexpr FormFactor H               = FormFactor(               constants::form_factor::H::a,               constants::form_factor::H::b,               constants::form_factor::H::c);
        constexpr FormFactor C               = FormFactor(               constants::form_factor::C::a,               constants::form_factor::C::b,               constants::form_factor::C::c);
        constexpr FormFactor N               = FormFactor(               constants::form_factor::N::a,               constants::form_factor::N::b,               constants::form_factor::N::c);
        constexpr FormFactor O               = FormFactor(               constants::form_factor::O::a,               constants::form_factor::O::b,               constants::form_factor::O::c);
        constexpr FormFactor S               = FormFactor(               constants::form_factor::S::a,               constants::form_factor::S::b,               constants::form_factor::S::c);

        // atomic groups
        constexpr FormFactor CH_sp3          = FormFactor(          constants::form_factor::CH_sp3::a,          constants::form_factor::CH_sp3::b,          constants::form_factor::CH_sp3::c);
        constexpr FormFactor CH2_sp3         = FormFactor(         constants::form_factor::CH2_sp3::a,         constants::form_factor::CH2_sp3::b,         constants::form_factor::CH2_sp3::c);
        constexpr FormFactor CH3_sp3         = FormFactor(         constants::form_factor::CH3_sp3::a,         constants::form_factor::CH3_sp3::b,         constants::form_factor::CH3_sp3::c);
        constexpr FormFactor CH_sp2          = FormFactor(          constants::form_factor::CH_sp2::a,          constants::form_factor::CH_sp2::b,          constants::form_factor::CH_sp2::c);
        constexpr FormFactor CH_arom         = FormFactor(         constants::form_factor::CH_arom::a,         constants::form_factor::CH_arom::b,         constants::form_factor::CH_arom::c);
        constexpr FormFactor OH_alc          = FormFactor(          constants::form_factor::OH_alc::a,          constants::form_factor::OH_alc::b,          constants::form_factor::OH_alc::c);
        constexpr FormFactor OH_acid         = FormFactor(         constants::form_factor::OH_acid::a,         constants::form_factor::OH_acid::b,         constants::form_factor::OH_acid::c);
        constexpr FormFactor O_res           = FormFactor(           constants::form_factor::O_res::a,           constants::form_factor::O_res::b,           constants::form_factor::O_res::c);
        constexpr FormFactor NH              = FormFactor(              constants::form_factor::NH::a,              constants::form_factor::NH::b,              constants::form_factor::NH::c);
        constexpr FormFactor NH2             = FormFactor(             constants::form_factor::NH2::a,             constants::form_factor::NH2::b,             constants::form_factor::NH2::c);
        constexpr FormFactor NH_plus         = FormFactor(         constants::form_factor::NH_plus::a,         constants::form_factor::NH_plus::b,         constants::form_factor::NH_plus::c);
        constexpr FormFactor NH2_plus        = FormFactor(        constants::form_factor::NH2_plus::a,        constants::form_factor::NH2_plus::b,        constants::form_factor::NH2_plus::c);
        constexpr FormFactor NH3_plus        = FormFactor(        constants::form_factor::NH3_plus::a,        constants::form_factor::NH3_plus::b,        constants::form_factor::NH3_plus::c);
        constexpr FormFactor NH_guanine      = FormFactor(      constants::form_factor::NH_guanine::a,      constants::form_factor::NH_guanine::b,      constants::form_factor::NH_guanine::c);
        constexpr FormFactor NH2_guanine     = FormFactor(     constants::form_factor::NH2_guanine::a,     constants::form_factor::NH2_guanine::b,     constants::form_factor::NH2_guanine::c);
        constexpr FormFactor SH              = FormFactor(              constants::form_factor::SH::a,              constants::form_factor::SH::b,              constants::form_factor::SH::c);

        // average excluded volume
        constexpr FormFactor excluded_volume = FormFactor( constants::form_factor::excluded_volume::a, constants::form_factor::excluded_volume::b, constants::form_factor::excluded_volume::c);
        // constexpr FormFactor excluded_volume = FormFactor(constants::displaced_volume::avg_vol);

        // all others; this is just the form factor of argon
        constexpr FormFactor other           = FormFactor(           constants::form_factor::other::a,           constants::form_factor::other::b,           constants::form_factor::other::c);

        constexpr const FormFactor& get_form_factor(form_factor_t type) {
            switch (type) {
                case form_factor_t::H:
                    return H;
                case form_factor_t::C:
                    return C;
                case form_factor_t::N:
                    return N;
                case form_factor_t::O:
                    return O;
                case form_factor_t::S:
                    return S;
                case form_factor_t::CH:
                    return CH_sp3;
                case form_factor_t::CH2:
                    return CH2_sp3;
                case form_factor_t::CH3:
                    return CH3_sp3;
                case form_factor_t::NH:
                    return NH;
                case form_factor_t::NH2:
                    return NH2;
                case form_factor_t::NH3:
                    return NH3_plus;
                case form_factor_t::OH:
                    return OH_alc;
                case form_factor_t::SH:
                    return SH;
                case form_factor_t::OTHER:
                    return other;
                case form_factor_t::EXCLUDED_VOLUME:
                    return excluded_volume;
                default:
                    throw std::runtime_error("form_factor::storage::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
        }
    };
}