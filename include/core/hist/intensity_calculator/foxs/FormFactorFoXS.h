// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <container/ArrayContainer2D.h>
#include <settings/FormFactorSettings.h>

#include <cmath>

namespace ausaxs::form_factor::foxs {
    class FormFactorFoXS {
        public:
            FormFactorFoXS(double q0) : q0(q0) {}

            double evaluate(double q) const {
                return q0*std::exp(-modulation*q*q);
            }

        private:
            double q0 = 1;
            double modulation = 0.23/2;
    };

    namespace storage {
        namespace atomic {
            static FormFactorFoXS get_form_factor(form_factor_t type) {
                switch (type) {
                    case form_factor_t::H:                  return FormFactorFoXS(0.999953);
                    case form_factor_t::C:                  return FormFactorFoXS(5.9992);
                    case form_factor_t::N:                  return FormFactorFoXS(6.9946);
                    case form_factor_t::O:                  return FormFactorFoXS(7.9994);
                    case form_factor_t::S:                  return FormFactorFoXS(15.9998);
                    case form_factor_t::CH:                 return FormFactorFoXS(6.99915);
                    case form_factor_t::CH2:                return FormFactorFoXS(7.99911);
                    case form_factor_t::CH3:                return FormFactorFoXS(8.99906);
                    case form_factor_t::NH:                 return FormFactorFoXS(7.99455);
                    case form_factor_t::NH2:                return FormFactorFoXS(8.99451);
                    case form_factor_t::NH3:                return FormFactorFoXS(9.99446);
                    case form_factor_t::OH:                 return FormFactorFoXS(8.99935);
                    case form_factor_t::SH:                 return FormFactorFoXS(16.9998);
                    case form_factor_t::OTHER:              return FormFactorFoXS(17.99);
                    case form_factor_t::EXCLUDED_VOLUME:    return FormFactorFoXS(0);
                    default:
                        throw std::runtime_error("form_factor::foxs::storage::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
            }

            [[maybe_unused]] static container::ArrayContainer2D<FormFactorProduct, settings::form_factor::max_ff_types, settings::form_factor::max_ff_types> generate_table() {
                auto ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
                container::ArrayContainer2D<FormFactorProduct, settings::form_factor::max_ff_types, settings::form_factor::max_ff_types> table;
                for (unsigned int i = form_factor::start_index_for_explicit_exv(); i < settings::form_factor::max_ff_types; ++i) {
                    for (unsigned int j = form_factor::start_index_for_explicit_exv(); j < i; ++j) {
                        table.index(i, j) = FormFactorProduct(
                            get_form_factor(static_cast<form_factor_t>(ff_indices[i])), 
                            get_form_factor(static_cast<form_factor_t>(ff_indices[j]))
                        );
                        table.index(j, i) = table.index(i, j);
                    }
                    table.index(i, i) = FormFactorProduct(
                        get_form_factor(static_cast<form_factor_t>(ff_indices[i])), 
                        get_form_factor(static_cast<form_factor_t>(ff_indices[i]))
                    );
                }
                return table;
            }
        }

        namespace exv {
            static FormFactorFoXS get_form_factor(form_factor_t type) {
                switch (type) {
                    case form_factor_t::H:      return FormFactorFoXS(1.7201);
                    case form_factor_t::C:      return FormFactorFoXS(5.49096);
                    case form_factor_t::N:      return FormFactorFoXS(0.83166);
                    case form_factor_t::O:      return FormFactorFoXS(3.04942);
                    case form_factor_t::S:      return FormFactorFoXS(6.63324);
                    case form_factor_t::CH:     return FormFactorFoXS(7.21106);
                    case form_factor_t::CH2:    return FormFactorFoXS(8.93116);
                    case form_factor_t::CH3:    return FormFactorFoXS(10.6513);
                    case form_factor_t::NH:     return FormFactorFoXS(2.55176);
                    case form_factor_t::NH2:    return FormFactorFoXS(4.27186);
                    case form_factor_t::NH3:    return FormFactorFoXS(5.99196);
                    case form_factor_t::OH:     return FormFactorFoXS(4.76952);
                    case form_factor_t::SH:     return FormFactorFoXS(8.35334);
                    case form_factor_t::OTHER:  return FormFactorFoXS(1.399);
                    default:
                        throw std::runtime_error("form_factor::foxs::storage::exv::get_form_factor: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
            }

            [[maybe_unused]] static container::ArrayContainer2D<FormFactorProduct, settings::form_factor::max_ff_types, settings::form_factor::max_ff_types> generate_table() {
                auto ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
                container::ArrayContainer2D<FormFactorProduct, settings::form_factor::max_ff_types, settings::form_factor::max_ff_types> table;
                for (unsigned int i = form_factor::start_index_for_explicit_exv(); i < settings::form_factor::max_ff_types; ++i) {
                    for (unsigned int j = form_factor::start_index_for_explicit_exv(); j < i; ++j) {
                        table.index(i, j) = FormFactorProduct(
                            get_form_factor(static_cast<form_factor_t>(ff_indices[i])), 
                            get_form_factor(static_cast<form_factor_t>(ff_indices[j]))
                        );
                        table.index(j, i) = table.index(i, j);
                    }
                    table.index(i, i) = FormFactorProduct(
                        get_form_factor(static_cast<form_factor_t>(ff_indices[i])), 
                        get_form_factor(static_cast<form_factor_t>(ff_indices[i]))
                    );
                }
                return table;
            }
        }

        namespace cross {
            [[maybe_unused]] static container::ArrayContainer2D<FormFactorProduct, settings::form_factor::max_ff_types, settings::form_factor::max_ff_types> generate_table() {
                auto ff_tables = form_factor::manager::get_active_product_tables();
                auto ff_indices = ff_tables->ff_indices;
                container::ArrayContainer2D<FormFactorProduct, settings::form_factor::max_ff_types, settings::form_factor::max_ff_types> table;
                for (unsigned int i = form_factor::start_index_for_explicit_exv(); i < ff_tables->active_count; ++i) {
                    for (unsigned int j = form_factor::start_index_for_explicit_exv(); j < ff_tables->active_count; ++j) {
                        table.index(i, j) = FormFactorProduct(
                            atomic::get_form_factor(static_cast<form_factor_t>(ff_indices[i])), 
                            exv::get_form_factor(static_cast<form_factor_t>(ff_indices[j]))
                        );
                    }
                }
                return table;
            }
        }
    }
}