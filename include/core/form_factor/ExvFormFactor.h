#pragma once

#include <constants/Constants.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/DisplacedVolumeTable.h>
#include <math/ConstexprMath.h>

namespace form_factor {
    /**
     * @brief Calculate the excluded volume form factor based on the description from Fraser, MacRae & Suzuki: https://doi.org/10.1107/S0021889878014296
     */
    class ExvFormFactor {
        public: 
            /**
             * @brief Create a new excluded volume form factor with the given volume.
             *
             * @param volume The excluded volume of the atom. 
             */
            constexpr ExvFormFactor(double volume) {
                exponent = constexpr_math::pow(volume, 2./3)/(4*constants::pi);
                q0 = volume*constants::charge::density::water;
            }

            constexpr double evaluate_normalized(double q) const {
                return constexpr_math::exp(-exponent*q*q);
            }

            constexpr double evaluate(double q) const {
                return q0*evaluate_normalized(q);
            }

            double exponent = 0;
            double q0 = 1;
    };

    namespace detail {
        struct ExvFormFactorSet {
            constexpr ExvFormFactorSet(const constants::displaced_volume::detail::DisplacedVolumeSet& set) : 
                H(set.H), C(set.C), CH(set.CH), CH2(set.CH2), CH3(set.CH3), 
                N(set.N), NH(set.NH), NH2(set.NH2), NH3(set.NH3), 
                O(set.O), OH(set.OH), 
                S(set.S), SH(set.SH),
                Ar(constants::displaced_volume::Ar)
            {}

            constexpr ExvFormFactor get_form_factor(form_factor_t type) const {
                switch(type) {
                    case form_factor_t::H:     return H;
                    case form_factor_t::C:     return C;
                    case form_factor_t::CH:    return CH;
                    case form_factor_t::CH2:   return CH2;
                    case form_factor_t::CH3:   return CH3;
                    case form_factor_t::N:     return N;
                    case form_factor_t::NH:    return NH;
                    case form_factor_t::NH2:   return NH2;
                    case form_factor_t::NH3:   return NH3;
                    case form_factor_t::O:     return O;
                    case form_factor_t::OH:    return OH;
                    case form_factor_t::S:     return S;
                    case form_factor_t::SH:    return SH;
                    case form_factor_t::OTHER: return Ar;
                    default: throw std::runtime_error("form_factor::detail::ExvFormFactorSet::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
            }

            ExvFormFactor H;
            ExvFormFactor C, CH, CH2, CH3;
            ExvFormFactor N, NH, NH2, NH3;
            ExvFormFactor O, OH;
            ExvFormFactor S, SH;
            ExvFormFactor Ar;
        };
    }

    namespace storage::exv {
        constexpr detail::ExvFormFactorSet standard(constants::displaced_volume::standard);
    }
}