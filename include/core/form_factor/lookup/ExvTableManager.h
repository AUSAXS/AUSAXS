// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/ExvTable.h>
#include <form_factor/FormFactorType.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::form_factor {
    class ExvTableManager {
        public:
            static const constants::exv::detail::ExvSet& get_default_exv_table();
            static observer_ptr<const constants::exv::detail::ExvSet> get_current_exv_table();

            /**
             * @brief Build the excluded volume form factors of the current set, displacing the solvent of the active probe.
             */
            static detail::ExvFormFactorSet get_current_exv_form_factor_set();

            static void set_custom_exv_table(const constants::exv::detail::ExvSet& set);
            static bool is_default();

            /**
             * @brief Get the total solvent volume displaced by @a molecule, in cubic angstroms.
             */
            static double get_total_displaced_volume(observer_ptr<const data::Molecule> molecule);

            /**
             * @brief Get the average solvent volume displaced by a single atom of @a molecule, in cubic angstroms.
             */
            static double get_average_displaced_volume(observer_ptr<const data::Molecule> molecule);

        private:
            static inline bool _use_custom_exv_table = false;
            static inline std::unique_ptr<constants::exv::detail::ExvSet> custom_exv_tables;
    };
}
