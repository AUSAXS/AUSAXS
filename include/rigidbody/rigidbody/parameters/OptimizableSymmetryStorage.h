// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/SymmetryStorage.h>

namespace ausaxs::symmetry {
    /**
     * @brief Marks a body's symmetries as owned outright, and thus safe for the optimiser to write to.
     *
     * A body holding ReferenceSymmetryViews shares its parameters with another body and must not be
     * converted, since writing through a view would silently update the owner's symmetry instead.
     */
    struct OptimizableSymmetryStorage : SymmetryStorage {
        OptimizableSymmetryStorage(SymmetryStorage&& other);
        ~OptimizableSymmetryStorage() override;

        std::unique_ptr<SymmetryStorage> clone() override;
    };
}
