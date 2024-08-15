#pragma once

#include <constants/Constants.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/DisplacedVolumeTable.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <container/ArrayContainer2D.h>
#include <math/ConstexprMath.h>

namespace form_factor::crysol {
    /**
     * @brief Calculate the excluded volume form factor based on the description from Crysol: https://doi.org/10.1107/S0021889895007047
     */
    struct ExvFormFactorCrysol {
        /**
            * @brief Create a new excluded volume form factor with the given volume.
            *
            * @param volume The excluded volume of the atom in cubic angstroms.
            */
        constexpr ExvFormFactorCrysol(double volume) {
            double magic_constant = 1/(4*constants::pi*constants::pi);
            exponent = magic_constant*constants::pi*constexpr_math::pow(volume, 2./3);
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

    namespace storage {
        struct exv {
            inline static ExvFormFactorCrysol H  =  ExvFormFactorCrysol(constants::displaced_volume::Traube.H);
            inline static ExvFormFactorCrysol C  =  ExvFormFactorCrysol(constants::displaced_volume::Traube.C);
            inline static ExvFormFactorCrysol CH  = ExvFormFactorCrysol(constants::displaced_volume::Traube.CH);
            inline static ExvFormFactorCrysol CH2 = ExvFormFactorCrysol(constants::displaced_volume::Traube.CH2);
            inline static ExvFormFactorCrysol CH3 = ExvFormFactorCrysol(constants::displaced_volume::Traube.CH3);
            inline static ExvFormFactorCrysol N  =  ExvFormFactorCrysol(constants::displaced_volume::Traube.N);
            inline static ExvFormFactorCrysol NH  = ExvFormFactorCrysol(constants::displaced_volume::Traube.NH);
            inline static ExvFormFactorCrysol NH2 = ExvFormFactorCrysol(constants::displaced_volume::Traube.NH2);
            inline static ExvFormFactorCrysol NH3 = ExvFormFactorCrysol(constants::displaced_volume::Traube.NH3);
            inline static ExvFormFactorCrysol O  =  ExvFormFactorCrysol(constants::displaced_volume::Traube.O);
            inline static ExvFormFactorCrysol OH  = ExvFormFactorCrysol(constants::displaced_volume::Traube.OH);
            inline static ExvFormFactorCrysol S  =  ExvFormFactorCrysol(constants::displaced_volume::Traube.S);
            inline static ExvFormFactorCrysol SH  = ExvFormFactorCrysol(constants::displaced_volume::Traube.SH);
            inline static ExvFormFactorCrysol Ar =  ExvFormFactorCrysol(constants::displaced_volume::Ar);

            inline static ExvFormFactorCrysol get_form_factor(form_factor_t type) {
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
                    default: throw std::runtime_error("form_factor::crysol::storage::exv::get_exv_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
            }

            [[maybe_unused]] static form_factor::storage::exv::table_t generate_table() {
                container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> table;
                for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
                    for (unsigned int j = 0; j < i; ++j) {
                        table.index(i, j) = PrecalculatedFormFactorProduct(
                            get_form_factor(static_cast<form_factor_t>(i)), 
                            get_form_factor(static_cast<form_factor_t>(j))
                        );
                        table.index(j, i) = table.index(i, j);
                    }
                    table.index(i, i) = PrecalculatedFormFactorProduct(
                        get_form_factor(static_cast<form_factor_t>(i)), 
                        get_form_factor(static_cast<form_factor_t>(i))
                    );
                }
                return table;
            }
        };

        struct cross {
            [[maybe_unused]] static form_factor::storage::cross::table_t generate_table() {
                container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> table;
                for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
                    for (unsigned int j = 0; j < i; ++j) {
                        table.index(i, j) = PrecalculatedFormFactorProduct(
                            form_factor::storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                            exv::get_form_factor(static_cast<form_factor_t>(j))
                        );
                        table.index(j, i) = table.index(i, j);
                    }
                    table.index(i, i) = PrecalculatedFormFactorProduct(
                        form_factor::storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                        exv::get_form_factor(static_cast<form_factor_t>(i))
                    );
                }
                return table;
            }
        };
    }
}