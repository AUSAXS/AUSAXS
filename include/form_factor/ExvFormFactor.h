#pragma once

#include <constants/Constants.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/DisplacedVolumeTable.h>

#include <cmath>

namespace form_factor {
    /**
     * @brief Calculate the excluded volume form factor based on the description from Fraser, MacRae & Suzuki: https://doi.org/10.1107/S0021889878014296
     *        This is the same approach as CRYSOL.
     */
    class ExvFormFactor {
        public: 
            constexpr ExvFormFactor(double volume) {
                q0 = volume*constants::charge::density::water;
                double s_to_q = constants::form_factor::s_to_q_factor; // crysol seems to have missed this conversion?
                exponent = M_PI*std::pow(volume, 2./3)*s_to_q;         // eq 6
            }

            constexpr double evaluate_normalized(double q) const {
                return std::exp(-exponent*q*q);
            }

            constexpr double evaluate(double q) const {
                return q0*evaluate_normalized(q);
            }

        private:
            double exponent = 0;
            double q0 = 1;
    };

    namespace storage::exv {
        constexpr ExvFormFactor CH =  ExvFormFactor(constants::DisplacedVolume::CH);
        constexpr ExvFormFactor CH2 = ExvFormFactor(constants::DisplacedVolume::CH2);
        constexpr ExvFormFactor CH3 = ExvFormFactor(constants::DisplacedVolume::CH3);
        constexpr ExvFormFactor NH =  ExvFormFactor(constants::DisplacedVolume::NH);
        constexpr ExvFormFactor NH2 = ExvFormFactor(constants::DisplacedVolume::NH2);
        constexpr ExvFormFactor NH3 = ExvFormFactor(constants::DisplacedVolume::NH3);
        constexpr ExvFormFactor OH =  ExvFormFactor(constants::DisplacedVolume::OH);
        constexpr ExvFormFactor SH =  ExvFormFactor(constants::DisplacedVolume::SH);

        constexpr ExvFormFactor H =  ExvFormFactor(constants::DisplacedVolume::H);
        constexpr ExvFormFactor C =  ExvFormFactor(constants::DisplacedVolume::C);
        constexpr ExvFormFactor N =  ExvFormFactor(constants::DisplacedVolume::N);
        constexpr ExvFormFactor O =  ExvFormFactor(constants::DisplacedVolume::O);
        constexpr ExvFormFactor S =  ExvFormFactor(constants::DisplacedVolume::S);
        constexpr ExvFormFactor Ar = ExvFormFactor(constants::DisplacedVolume::Ar);

        constexpr ExvFormFactor get_form_factor(form_factor_t type) {
            switch(type) {
                case form_factor_t::CH: return CH;
                case form_factor_t::CH2: return CH2;
                case form_factor_t::CH3: return CH3;
                case form_factor_t::NH: return NH;
                case form_factor_t::NH2: return NH2;
                case form_factor_t::OH: return OH;
                case form_factor_t::H: return H;
                case form_factor_t::C: return C;
                case form_factor_t::N: return N;
                case form_factor_t::O: return O;
                case form_factor_t::S: return S;
                case form_factor_t::SH: return SH;
                case form_factor_t::OTHER: return Ar;
                default: throw std::runtime_error("form_factor::storage::exv::get_exv_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
        }
    }
}