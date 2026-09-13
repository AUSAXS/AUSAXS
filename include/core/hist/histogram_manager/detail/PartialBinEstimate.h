// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>

namespace ausaxs::hist::detail {
    /**
     * @brief Return a conservative bin count for a persistent partial-manager axis.
     *
     * The bound includes the molecule's atoms and waters, plus transformed corners of
     * each symmetric copy that has not been materialized in the molecule.
     */
    template<bool variable_bin_width>
    int required_partial_bin_count(const data::Molecule& protein);

    /**
     * @brief Add capacity for ordinary movement between partial-manager rebuilds.
     */
    int grown_partial_bin_count(int required);
}