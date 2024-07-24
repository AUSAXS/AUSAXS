#pragma once

#include <constants/Constants.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/DisplacedVolumeTable.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <container/ArrayContainer2D.h>
#include <math/ConstexprMath.h>

#include <algorithm>

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

    struct FormFactorCrysol {
        constexpr FormFactorCrysol(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {
            std::transform(b.begin(), b.end(), this->b.begin(), [](double x) { return x; });
            f0 = 1./(a[0] + a[1] + a[2] + a[3] + a[4] + c);
        }

        constexpr double evaluate(double q) const {
            double sum = 0;
            for (unsigned int i = 0; i < 5; ++i) {
                sum += a[i]*constexpr_math::exp(-b[i]*q*q);
            }
            return (sum + c)*f0;
        }

        std::array<double, 5> a;
        std::array<double, 5> b;
        double c;
        double f0 = 1;
    };

    namespace storage {
        struct atomic {
            inline static FormFactorCrysol H            = FormFactorCrysol(           constants::form_factor::H::a,           constants::form_factor::H::b,           constants::form_factor::H::c);
            inline static FormFactorCrysol C            = FormFactorCrysol(           constants::form_factor::C::a,           constants::form_factor::C::b,           constants::form_factor::C::c);
            inline static FormFactorCrysol N            = FormFactorCrysol(           constants::form_factor::N::a,           constants::form_factor::N::b,           constants::form_factor::N::c);
            inline static FormFactorCrysol O            = FormFactorCrysol(           constants::form_factor::O::a,           constants::form_factor::O::b,           constants::form_factor::O::c);
            inline static FormFactorCrysol S            = FormFactorCrysol(           constants::form_factor::S::a,           constants::form_factor::S::b,           constants::form_factor::S::c);

            inline static FormFactorCrysol CH_sp3       = FormFactorCrysol(      constants::form_factor::CH_sp3::a,      constants::form_factor::CH_sp3::b,      constants::form_factor::CH_sp3::c);
            inline static FormFactorCrysol CH2_sp3      = FormFactorCrysol(     constants::form_factor::CH2_sp3::a,     constants::form_factor::CH2_sp3::b,     constants::form_factor::CH2_sp3::c);
            inline static FormFactorCrysol CH3_sp3      = FormFactorCrysol(     constants::form_factor::CH3_sp3::a,     constants::form_factor::CH3_sp3::b,     constants::form_factor::CH3_sp3::c);
            inline static FormFactorCrysol CH_sp2       = FormFactorCrysol(      constants::form_factor::CH_sp2::a,      constants::form_factor::CH_sp2::b,      constants::form_factor::CH_sp2::c);
            inline static FormFactorCrysol CH_arom      = FormFactorCrysol(     constants::form_factor::CH_arom::a,     constants::form_factor::CH_arom::b,     constants::form_factor::CH_arom::c);
            inline static FormFactorCrysol OH_alc       = FormFactorCrysol(      constants::form_factor::OH_alc::a,      constants::form_factor::OH_alc::b,      constants::form_factor::OH_alc::c);
            inline static FormFactorCrysol OH_acid      = FormFactorCrysol(     constants::form_factor::OH_acid::a,     constants::form_factor::OH_acid::b,     constants::form_factor::OH_acid::c);
            inline static FormFactorCrysol O_res        = FormFactorCrysol(       constants::form_factor::O_res::a,       constants::form_factor::O_res::b,       constants::form_factor::O_res::c);
            inline static FormFactorCrysol NH           = FormFactorCrysol(          constants::form_factor::NH::a,          constants::form_factor::NH::b,          constants::form_factor::NH::c);
            inline static FormFactorCrysol NH2          = FormFactorCrysol(         constants::form_factor::NH2::a,         constants::form_factor::NH2::b,         constants::form_factor::NH2::c);
            inline static FormFactorCrysol NH_plus      = FormFactorCrysol(     constants::form_factor::NH_plus::a,     constants::form_factor::NH_plus::b,     constants::form_factor::NH_plus::c);
            inline static FormFactorCrysol NH2_plus     = FormFactorCrysol(    constants::form_factor::NH2_plus::a,    constants::form_factor::NH2_plus::b,    constants::form_factor::NH2_plus::c);
            inline static FormFactorCrysol NH3_plus     = FormFactorCrysol(    constants::form_factor::NH3_plus::a,    constants::form_factor::NH3_plus::b,    constants::form_factor::NH3_plus::c);
            inline static FormFactorCrysol NH_guanine   = FormFactorCrysol(  constants::form_factor::NH_guanine::a,  constants::form_factor::NH_guanine::b,  constants::form_factor::NH_guanine::c);
            inline static FormFactorCrysol NH2_guanine  = FormFactorCrysol( constants::form_factor::NH2_guanine::a, constants::form_factor::NH2_guanine::b, constants::form_factor::NH2_guanine::c);
            inline static FormFactorCrysol SH           = FormFactorCrysol(          constants::form_factor::SH::a,          constants::form_factor::SH::b,          constants::form_factor::SH::c);

            inline static FormFactorCrysol get_form_factor(form_factor_t type) {
                switch(type) {
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
                    case form_factor_t::OTHER:              return O_res;
                    case form_factor_t::EXCLUDED_VOLUME:    return FormFactorCrysol({0,0,0,0,0}, {0,0,0,0,0}, 0); // unused
                    default: throw std::runtime_error("form_factor::crysol::storage::atomic::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
            }

            [[maybe_unused]] static form_factor::storage::atomic::table_t generate_table() {
                container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()> table;
                for (unsigned int i = 0; i < form_factor::get_count(); ++i) {
                    for (unsigned int j = 0; j < i; ++j) {
                        table.index(i, j) = PrecalculatedFormFactorProduct(
                            form_factor::crysol::storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                            form_factor::crysol::storage::atomic::get_form_factor(static_cast<form_factor_t>(j))
                        );
                        table.index(j, i) = table.index(i, j);
                    }
                    table.index(i, i) = PrecalculatedFormFactorProduct(
                        form_factor::crysol::storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                        form_factor::crysol::storage::atomic::get_form_factor(static_cast<form_factor_t>(i))
                    );
                }
                return table;
            }
        };

        struct exv {
            inline static ExvFormFactorCrysol CH  = ExvFormFactorCrysol(constants::displaced_volume::CH);
            inline static ExvFormFactorCrysol CH2 = ExvFormFactorCrysol(constants::displaced_volume::CH2);
            inline static ExvFormFactorCrysol CH3 = ExvFormFactorCrysol(constants::displaced_volume::CH3);
            inline static ExvFormFactorCrysol NH  = ExvFormFactorCrysol(constants::displaced_volume::NH);
            inline static ExvFormFactorCrysol NH2 = ExvFormFactorCrysol(constants::displaced_volume::NH2);
            inline static ExvFormFactorCrysol NH3 = ExvFormFactorCrysol(constants::displaced_volume::NH3);
            inline static ExvFormFactorCrysol OH  = ExvFormFactorCrysol(constants::displaced_volume::OH);
            inline static ExvFormFactorCrysol SH  = ExvFormFactorCrysol(constants::displaced_volume::SH);

            inline static ExvFormFactorCrysol H  =  ExvFormFactorCrysol(constants::displaced_volume::H);
            inline static ExvFormFactorCrysol C  =  ExvFormFactorCrysol(constants::displaced_volume::C);
            inline static ExvFormFactorCrysol N  =  ExvFormFactorCrysol(constants::displaced_volume::N);
            inline static ExvFormFactorCrysol O  =  ExvFormFactorCrysol(constants::displaced_volume::O);
            inline static ExvFormFactorCrysol S  =  ExvFormFactorCrysol(constants::displaced_volume::S);
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
                            atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                            exv::get_form_factor(static_cast<form_factor_t>(j))
                        );
                        table.index(j, i) = table.index(i, j);
                    }
                    table.index(i, i) = PrecalculatedFormFactorProduct(
                        atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                        exv::get_form_factor(static_cast<form_factor_t>(i))
                    );
                }
                return table;
            }
        };
    }
}