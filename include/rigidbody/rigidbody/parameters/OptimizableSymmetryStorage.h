// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/SymmetryStorage.h>

namespace ausaxs::symmetry {
    struct OptimizableSymmetryStorage : SymmetryStorage {
        OptimizableSymmetryStorage(SymmetryStorage&& other);
        ~OptimizableSymmetryStorage() override;

        void add(symmetry::type symmetry) override;

        std::unique_ptr<SymmetryStorage> clone() override;

        bool optimize_translate = false;
        bool optimize_rot_axis = false;
    };
}