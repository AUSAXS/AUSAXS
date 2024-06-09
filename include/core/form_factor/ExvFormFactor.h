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

    namespace storage::exv {
        constexpr ExvFormFactor CH  = ExvFormFactor(constants::displaced_volume::CH);
        constexpr ExvFormFactor CH2 = ExvFormFactor(constants::displaced_volume::CH2);
        constexpr ExvFormFactor CH3 = ExvFormFactor(constants::displaced_volume::CH3);
        constexpr ExvFormFactor NH  = ExvFormFactor(constants::displaced_volume::NH);
        constexpr ExvFormFactor NH2 = ExvFormFactor(constants::displaced_volume::NH2);
        constexpr ExvFormFactor NH3 = ExvFormFactor(constants::displaced_volume::NH3);
        constexpr ExvFormFactor OH  = ExvFormFactor(constants::displaced_volume::OH);
        constexpr ExvFormFactor SH  = ExvFormFactor(constants::displaced_volume::SH);

        constexpr ExvFormFactor H  =  ExvFormFactor(constants::displaced_volume::H);
        constexpr ExvFormFactor C  =  ExvFormFactor(constants::displaced_volume::C);
        constexpr ExvFormFactor N  =  ExvFormFactor(constants::displaced_volume::N);
        constexpr ExvFormFactor O  =  ExvFormFactor(constants::displaced_volume::O);
        constexpr ExvFormFactor S  =  ExvFormFactor(constants::displaced_volume::S);
        constexpr ExvFormFactor Ar =  ExvFormFactor(constants::displaced_volume::Ar);

        constexpr ExvFormFactor get_form_factor(form_factor_t type) {
            switch(type) {
                case form_factor_t::CH:    return CH;
                case form_factor_t::CH2:   return CH2;
                case form_factor_t::CH3:   return CH3;
                case form_factor_t::NH:    return NH;
                case form_factor_t::NH2:   return NH2;
                case form_factor_t::NH3:   return NH3;
                case form_factor_t::OH:    return OH;
                case form_factor_t::H:     return H;
                case form_factor_t::C:     return C;
                case form_factor_t::N:     return N;
                case form_factor_t::O:     return O;
                case form_factor_t::S:     return S;
                case form_factor_t::SH:    return SH;
                case form_factor_t::OTHER: return Ar;
                default: throw std::runtime_error("form_factor::storage::exv::get_exv_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
        }
    }
}