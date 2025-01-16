#pragma once

#include <data/symmetry/SymmetryStorage.h>

namespace ausaxs::symmetry  {
    struct OptimizableSymmetryStorage : SymmetryStorage {
        virtual ~OptimizableSymmetryStorage() = default;

        virtual std::unique_ptr<SymmetryStorage> clone();

        bool optimize_translate = false;
        bool optimize_rotate_cm = false;
    };
}

inline std::unique_ptr<ausaxs::symmetry::SymmetryStorage> ausaxs::symmetry::OptimizableSymmetryStorage::clone() {
    return std::make_unique<ausaxs::symmetry::OptimizableSymmetryStorage>(*static_cast<ausaxs::symmetry::OptimizableSymmetryStorage*>(this));
}