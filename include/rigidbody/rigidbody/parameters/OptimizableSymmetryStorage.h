#pragma once

#include <data/symmetry/SymmetryStorage.h>

namespace ausaxs::symmetry  {
    struct OptimizableSymmetryStorage : SymmetryStorage {
        ~OptimizableSymmetryStorage() override = default;

        void add(symmetry::type symmetry) override;

        std::unique_ptr<SymmetryStorage> clone() override;

        bool optimize_translate = false;
        bool optimize_rotate_cm = false;
    };
}

inline std::unique_ptr<ausaxs::symmetry::SymmetryStorage> ausaxs::symmetry::OptimizableSymmetryStorage::clone() {
    return std::make_unique<ausaxs::symmetry::OptimizableSymmetryStorage>(*static_cast<ausaxs::symmetry::OptimizableSymmetryStorage*>(this));
}

inline void ausaxs::symmetry::OptimizableSymmetryStorage::add(symmetry::type symmetry) {
    symmetries.emplace_back(symmetry::get(symmetry));
    switch (symmetry) {
        case symmetry::type::p2:
        case symmetry::type::p3:
        case symmetry::type::p4:
            optimize_rotate_cm = true;
            break;

        default:
            throw std::runtime_error("OptimizableSymmetryStorage::add: Symmetry \"" + std::to_string(static_cast<int>(symmetry)) + "\" not found. Did you forget to add it?");
    }
}