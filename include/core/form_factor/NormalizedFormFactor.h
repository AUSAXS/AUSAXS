// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactor.h>

#include <stdexcept>

namespace ausaxs::form_factor {
    struct NormalizedFormFactor : public FormFactor {
        constexpr NormalizedFormFactor(ExvFormFactor&& ffx) : FormFactor(std::move(ffx)) {set_normalization(1);}
        constexpr NormalizedFormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : FormFactor(a, b, c) {set_normalization(1);}
        constexpr NormalizedFormFactor(FormFactor&& ff) : FormFactor(std::move(ff)) {set_normalization(1);}
        constexpr NormalizedFormFactor(const FormFactor& ff) : FormFactor(ff) {set_normalization(1);}
    };

    /**
     * This struct contains the form factors of the most common atomic elements encountered in SAXS. 
     */
    namespace lookup::atomic::normalized {
        // atomic
        constexpr NormalizedFormFactor H               = NormalizedFormFactor(               constants::form_factor::H::a,               constants::form_factor::H::b,               constants::form_factor::H::c);
        constexpr NormalizedFormFactor C               = NormalizedFormFactor(               constants::form_factor::C::a,               constants::form_factor::C::b,               constants::form_factor::C::c);
        constexpr NormalizedFormFactor N               = NormalizedFormFactor(               constants::form_factor::N::a,               constants::form_factor::N::b,               constants::form_factor::N::c);
        constexpr NormalizedFormFactor O               = NormalizedFormFactor(               constants::form_factor::O::a,               constants::form_factor::O::b,               constants::form_factor::O::c);
        constexpr NormalizedFormFactor S               = NormalizedFormFactor(               constants::form_factor::S::a,               constants::form_factor::S::b,               constants::form_factor::S::c);

        // atomic groups
        constexpr NormalizedFormFactor CH_sp3          = NormalizedFormFactor(          constants::form_factor::CH_sp3::a,          constants::form_factor::CH_sp3::b,          constants::form_factor::CH_sp3::c);
        constexpr NormalizedFormFactor CH2_sp3         = NormalizedFormFactor(         constants::form_factor::CH2_sp3::a,         constants::form_factor::CH2_sp3::b,         constants::form_factor::CH2_sp3::c);
        constexpr NormalizedFormFactor CH3_sp3         = NormalizedFormFactor(         constants::form_factor::CH3_sp3::a,         constants::form_factor::CH3_sp3::b,         constants::form_factor::CH3_sp3::c);
        constexpr NormalizedFormFactor CH_sp2          = NormalizedFormFactor(          constants::form_factor::CH_sp2::a,          constants::form_factor::CH_sp2::b,          constants::form_factor::CH_sp2::c);
        constexpr NormalizedFormFactor CH_arom         = NormalizedFormFactor(         constants::form_factor::CH_arom::a,         constants::form_factor::CH_arom::b,         constants::form_factor::CH_arom::c);
        constexpr NormalizedFormFactor OH_alc          = NormalizedFormFactor(          constants::form_factor::OH_alc::a,          constants::form_factor::OH_alc::b,          constants::form_factor::OH_alc::c);
        constexpr NormalizedFormFactor OH_acid         = NormalizedFormFactor(         constants::form_factor::OH_acid::a,         constants::form_factor::OH_acid::b,         constants::form_factor::OH_acid::c);
        constexpr NormalizedFormFactor O_res           = NormalizedFormFactor(           constants::form_factor::O_res::a,           constants::form_factor::O_res::b,           constants::form_factor::O_res::c);
        constexpr NormalizedFormFactor NH              = NormalizedFormFactor(              constants::form_factor::NH::a,              constants::form_factor::NH::b,              constants::form_factor::NH::c);
        constexpr NormalizedFormFactor NH2             = NormalizedFormFactor(             constants::form_factor::NH2::a,             constants::form_factor::NH2::b,             constants::form_factor::NH2::c);
        constexpr NormalizedFormFactor NH_plus         = NormalizedFormFactor(         constants::form_factor::NH_plus::a,         constants::form_factor::NH_plus::b,         constants::form_factor::NH_plus::c);
        constexpr NormalizedFormFactor NH2_plus        = NormalizedFormFactor(        constants::form_factor::NH2_plus::a,        constants::form_factor::NH2_plus::b,        constants::form_factor::NH2_plus::c);
        constexpr NormalizedFormFactor NH3_plus        = NormalizedFormFactor(        constants::form_factor::NH3_plus::a,        constants::form_factor::NH3_plus::b,        constants::form_factor::NH3_plus::c);
        constexpr NormalizedFormFactor NH_guanine      = NormalizedFormFactor(      constants::form_factor::NH_guanine::a,      constants::form_factor::NH_guanine::b,      constants::form_factor::NH_guanine::c);
        constexpr NormalizedFormFactor NH2_guanine     = NormalizedFormFactor(     constants::form_factor::NH2_guanine::a,     constants::form_factor::NH2_guanine::b,     constants::form_factor::NH2_guanine::c);
        constexpr NormalizedFormFactor SH              = NormalizedFormFactor(              constants::form_factor::SH::a,              constants::form_factor::SH::b,              constants::form_factor::SH::c);

        // average excluded volume
        constexpr NormalizedFormFactor excluded_volume = NormalizedFormFactor( constants::form_factor::excluded_volume::a, constants::form_factor::excluded_volume::b, constants::form_factor::excluded_volume::c);

        // all others; this is just the form factor of argon
        constexpr NormalizedFormFactor other           = NormalizedFormFactor(           constants::form_factor::other::a,           constants::form_factor::other::b,           constants::form_factor::other::c);

        constexpr const NormalizedFormFactor& get(form_factor_t type) {
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
                case form_factor_t::UNKNOWN:
                throw std::runtime_error(
                    "form_factor::lookup::atomic::normalized::get: Attempted to get the form factor of an UNKNOWN atom.\n"
                    "This typically occurs when performing species-dependent operations on data without form factor information."
                );
                default: throw std::runtime_error("form_factor::lookup::atomic::normalized::get: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
        }
    }
}
