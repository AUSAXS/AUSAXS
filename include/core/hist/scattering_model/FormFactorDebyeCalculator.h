// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>

#include <vector>

namespace ausaxs::hist {
    /**
     * @brief The Debye intensity of the atoms and waters of a molecule with their atomic form factors, without any excluded volume:
     *      I(q) = \sum_{i != j} f_i(q) f_j(q) sinc(q r_ij) + \sum_i s_i(q),
     *        where f is the amplitude of the form factor type of each scatterer, and s its self-correlation form factor.
     *        The radiation is that of settings::scattering::radiation.
     *
     * Every form factor type present in the molecule is resolved, irrespective of settings::form_factor, so rare elements such as a
     * single metal ion keep their own form factor. The form factors are evaluated at each q value, so any q is supported, also
     * beyond the default q axis. The distances are binned with weighted bins.
     *
     * The weights of the atoms are not used; each scatterer contributes the form factor of its type.
     * Throws if any atom lacks form factor information, or if the molecule has symmetries.
     */
    std::vector<double> form_factor_debye_transform(const data::Molecule& molecule, const std::vector<double>& q);
}
