#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <container/ArrayContainer2D.h>

#include <cmath>

namespace form_factor::foxs {
    class FormFactorFoXS {
        public:
            FormFactorFoXS(double q0) : q0(q0) {}
            FormFactorFoXS(double q0, double) : q0(q0) {}
            // FormFactorFoXS(double q0, double modulation) : q0(q0), modulation(modulation) {}

            double evaluate(double q) const {
                return q0*std::exp(-modulation*q*q);
            }

        private:
            double q0 = 1;
            double modulation = 0.23/2;
    };
    static_assert(sizeof(FormFactorFoXS) == sizeof(ExvFormFactor), "ExvFormFactorCrysol must have the same size as ExvFormFactor");

    namespace storage {
        namespace atomic {
            static double vacuum_modulation = 0.1/2; 
            static FormFactorFoXS get_form_factor(form_factor_t type) {
                switch (type) {
                    case form_factor_t::H:                  return FormFactorFoXS(0.999953, vacuum_modulation);
                    case form_factor_t::C:                  return FormFactorFoXS(5.9992,   vacuum_modulation);
                    case form_factor_t::N:                  return FormFactorFoXS(6.9946,   vacuum_modulation);
                    case form_factor_t::O:                  return FormFactorFoXS(7.9994,   vacuum_modulation);
                    case form_factor_t::S:                  return FormFactorFoXS(15.9998,  vacuum_modulation);
                    case form_factor_t::CH:                 return FormFactorFoXS(6.99915,  vacuum_modulation);
                    case form_factor_t::CH2:                return FormFactorFoXS(7.99911,  vacuum_modulation);
                    case form_factor_t::CH3:                return FormFactorFoXS(8.99906,  vacuum_modulation);
                    case form_factor_t::NH:                 return FormFactorFoXS(7.99455,  vacuum_modulation);
                    case form_factor_t::NH2:                return FormFactorFoXS(8.99451,  vacuum_modulation);
                    case form_factor_t::NH3:                return FormFactorFoXS(9.99446,  vacuum_modulation);
                    case form_factor_t::OH:                 return FormFactorFoXS(8.99935,  vacuum_modulation);
                    case form_factor_t::SH:                 return FormFactorFoXS(16.9998,  vacuum_modulation);
                    case form_factor_t::OTHER:              return FormFactorFoXS(17.99,    vacuum_modulation);
                    case form_factor_t::EXCLUDED_VOLUME:    return FormFactorFoXS(0,        vacuum_modulation);
                    default:
                        throw std::runtime_error("form_factor::foxs::storage::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
            }

            [[maybe_unused]] static container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()> generate_table() {
                container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()> table;
                for (unsigned int i = 0; i < form_factor::get_count(); ++i) {
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
        }

        namespace exv {
            static double exv_modulation = 1.62*1.62/2; 
            static FormFactorFoXS get_form_factor(form_factor_t type) {
                switch (type) {
                    case form_factor_t::H:      return FormFactorFoXS(1.7201,   exv_modulation);
                    case form_factor_t::C:      return FormFactorFoXS(5.49096,  exv_modulation);
                    case form_factor_t::N:      return FormFactorFoXS(0.83166,  exv_modulation);
                    case form_factor_t::O:      return FormFactorFoXS(3.04942,  exv_modulation);
                    case form_factor_t::S:      return FormFactorFoXS(6.63324,  exv_modulation);
                    case form_factor_t::CH:     return FormFactorFoXS(7.21106,  exv_modulation);
                    case form_factor_t::CH2:    return FormFactorFoXS(8.93116,  exv_modulation);
                    case form_factor_t::CH3:    return FormFactorFoXS(10.6513,  exv_modulation);
                    case form_factor_t::NH:     return FormFactorFoXS(2.55176,  exv_modulation);
                    case form_factor_t::NH2:    return FormFactorFoXS(4.27186,  exv_modulation);
                    case form_factor_t::NH3:    return FormFactorFoXS(5.99196,  exv_modulation);
                    case form_factor_t::OH:     return FormFactorFoXS(4.76952,  exv_modulation);
                    case form_factor_t::SH:     return FormFactorFoXS(8.35334,  exv_modulation);
                    case form_factor_t::OTHER:  return FormFactorFoXS(1.399,    exv_modulation);
                    default:
                        throw std::runtime_error("form_factor::foxs::storage::exv::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
            }

            [[maybe_unused]] static container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> generate_table() {
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
        }

        namespace cross {
            [[maybe_unused]] static container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> generate_table() {
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
        }
    }
}