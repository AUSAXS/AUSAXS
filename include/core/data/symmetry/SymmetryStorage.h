// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <data/symmetry/Symmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <vector>

namespace ausaxs::symmetry {
    struct SymmetryStorage {
        SymmetryStorage() = default;
        SymmetryStorage(const SymmetryStorage&) = delete;
        SymmetryStorage& operator=(const SymmetryStorage&) = delete;
        SymmetryStorage(SymmetryStorage&&) = default;
        SymmetryStorage& operator=(SymmetryStorage&&) = default;
        virtual ~SymmetryStorage() = default;

        observer_ptr<ISymmetry> get(std::size_t i);
        observer_ptr<const ISymmetry> get(std::size_t i) const;

        const std::vector<std::unique_ptr<ISymmetry>>& get() const;
        std::vector<std::unique_ptr<ISymmetry>>& get();

        virtual void add(symmetry::type symmetry);
        void add(std::unique_ptr<ISymmetry> sym);

        virtual std::unique_ptr<SymmetryStorage> clone();

        std::vector<std::unique_ptr<ISymmetry>> symmetries;
    };
}

inline ausaxs::observer_ptr<ausaxs::symmetry::ISymmetry> ausaxs::symmetry::SymmetryStorage::get(std::size_t i) {
    return symmetries[i].get();
}

inline ausaxs::observer_ptr<const ausaxs::symmetry::ISymmetry> ausaxs::symmetry::SymmetryStorage::get(std::size_t i) const {
    return symmetries[i].get();
}

inline const std::vector<std::unique_ptr<ausaxs::symmetry::ISymmetry>>& ausaxs::symmetry::SymmetryStorage::get() const {
    return symmetries;
}

inline std::vector<std::unique_ptr<ausaxs::symmetry::ISymmetry>>& ausaxs::symmetry::SymmetryStorage::get() {
    return symmetries;
}

inline std::unique_ptr<ausaxs::symmetry::SymmetryStorage> ausaxs::symmetry::SymmetryStorage::clone() {
    auto copy = std::make_unique<SymmetryStorage>();
    for (const auto& s : symmetries) { copy->symmetries.push_back(s->clone()); }
    return copy;
}

inline void ausaxs::symmetry::SymmetryStorage::add(symmetry::type symmetry) {
    symmetries.push_back(symmetry::get(symmetry));
}

inline void ausaxs::symmetry::SymmetryStorage::add(std::unique_ptr<ISymmetry> sym) {
    symmetries.push_back(std::move(sym));
}
