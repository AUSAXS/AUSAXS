// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorInfo.h>
#include <form_factor/GroupFormFactor.h>
#include <form_factor/neutron/NeutronFormFactorTable.h>

#include <cassert>

namespace ausaxs::form_factor::neutron {
    /**
     * @brief The neutron form factor of a single nucleus: its coherent scattering length, which does not depend on q.
     */
    class AtomicFormFactor {
        public:
            constexpr AtomicFormFactor(double scattering_length) : b(scattering_length) {}

            /**
             * @brief Evaluate the form factor at a given q value.
             */
            constexpr double evaluate(double /*q*/) const {return b;}

            /**
             * @brief Evaluate the form factor at q = 0.
             */
            constexpr double I0() const {return b;}

        private:
            double b;
    };

    /**
     * @brief The neutron form factor of a form factor type: a heavy atom with its bound hydrogens.
     */
    using FormFactor = GroupFormFactor<AtomicFormFactor>;

    /**
     * @brief The deuterium fraction of the implicit hydrogens of a form factor type.
     *
     * @param d2o_fraction The D2O fraction of the solvent. It may lie outside [0, 1], see neutron::d2o_fraction.
     * @param exchange_efficiency The fraction of the exchangeable hydrogens that has exchanged with the solvent.
     */
    constexpr double deuterium_fraction(form_factor_t type, double d2o_fraction, double exchange_efficiency) {
        if (type == form_factor_t::WATER) {return d2o_fraction;} // this is just the solvent composition
        switch (form_factor::detail::ff_info_table[static_cast<int>(type)].element) {
            case constants::atom_t::N:
            case constants::atom_t::O:
            case constants::atom_t::S:
                return exchange_efficiency*d2o_fraction;
            default: // C and other elements do not exchange, so their implicit hydrogens are always protium
                return 0;
        }
    }

    /**
     * @brief Get the neutron form factor of a form factor type.
     *        Each implicit hydrogen has the mean scattering length (1-x)*b_H + x*b_D, with x its deuterium_fraction.
     *
     * @param d2o_fraction The D2O fraction of the solvent.
     * @param exchange_efficiency The fraction of the exchangeable hydrogens that has exchanged with the solvent.
     */
    constexpr FormFactor get(form_factor_t type, double d2o_fraction, double exchange_efficiency) {
        assert(form_factor::detail::is_tabulated(type) && "form_factor::neutron::get: Invalid form factor type. Untabulated types must be rejected before reaching here.");
        assert(type != form_factor_t::EXCLUDED_VOLUME && "form_factor::neutron::get: The excluded volume form factor is not defined for neutrons.");
        const auto& info = form_factor::detail::ff_info_table[static_cast<int>(type)];
        const double x = deuterium_fraction(type, d2o_fraction, exchange_efficiency);
        return {
            AtomicFormFactor(info.neutron_scattering_length),
            AtomicFormFactor((1-x)*scattering_length::H + x*scattering_length::D),
            info.geometry
        };
    }

    /**
     * @brief Get the neutron form factor of a form factor type in the solvent as prescibed by the current settings.
     */
    FormFactor get(form_factor_t type);
}
