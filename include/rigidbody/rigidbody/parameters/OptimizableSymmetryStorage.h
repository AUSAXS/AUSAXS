// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/SymmetryStorage.h>

namespace ausaxs::symmetry  {
    struct OptimizableSymmetryStorage : SymmetryStorage {
        OptimizableSymmetryStorage(SymmetryStorage&& other) : SymmetryStorage(std::move(other)) {}
        ~OptimizableSymmetryStorage() override = default;

        void add(symmetry::type symmetry) override;

        std::unique_ptr<SymmetryStorage> clone() override;

        bool optimize_translate = false;
        bool optimize_rot_axis = false;
        bool optimize_rot_angle = false;
    };
}

inline std::unique_ptr<ausaxs::symmetry::SymmetryStorage> ausaxs::symmetry::OptimizableSymmetryStorage::clone() {
    return std::make_unique<ausaxs::symmetry::OptimizableSymmetryStorage>(*static_cast<ausaxs::symmetry::OptimizableSymmetryStorage*>(this));
}

inline void ausaxs::symmetry::OptimizableSymmetryStorage::add(symmetry::type symmetry) {
    symmetries.emplace_back(symmetry::get(symmetry));
    switch (symmetry) {
        case symmetry::type::c2:
        case symmetry::type::c3:
        case symmetry::type::c4:
        case symmetry::type::c5:
        case symmetry::type::c6:
            optimize_rot_axis = true;
            optimize_translate = true;
            break;

        case symmetry::type::p2:
            optimize_rot_axis = true;
            optimize_translate = true;
            optimize_rot_angle = true;
            break;

        default:
            throw std::runtime_error("OptimizableSymmetryStorage::add: Symmetry \"" + std::to_string(static_cast<int>(symmetry)) + "\" not found. Did you forget to add it?");
    }
}