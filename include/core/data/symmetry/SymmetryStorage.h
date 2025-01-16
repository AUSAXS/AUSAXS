#pragma once

#include <data/symmetry/Symmetry.h>

namespace ausaxs::symmetry {
    struct SymmetryStorage {
        virtual ~SymmetryStorage() = default;

        const std::vector<Symmetry>& get() const;
        std::vector<Symmetry>& get();

        virtual std::unique_ptr<SymmetryStorage> clone();

        std::vector<Symmetry> symmetries;
    };
}

inline const std::vector<ausaxs::symmetry::Symmetry>& ausaxs::symmetry::SymmetryStorage::get() const {
    return symmetries;
}

inline std::vector<ausaxs::symmetry::Symmetry>& ausaxs::symmetry::SymmetryStorage::get() {
    return symmetries;
}

inline std::unique_ptr<ausaxs::symmetry::SymmetryStorage> ausaxs::symmetry::SymmetryStorage::clone() {
    return std::make_unique<SymmetryStorage>(*this);
}