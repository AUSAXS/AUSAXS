// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <data/symmetry/ReferenceSymmetry.h>

#include <cassert>

using namespace ausaxs::symmetry;

OptimizableSymmetryStorage::OptimizableSymmetryStorage(SymmetryStorage&& other) : SymmetryStorage(std::move(other)) {
    assert([this] () {
        for (const auto& s : this->symmetries) {
            if (dynamic_cast<const ReferenceSymmetryView*>(s.get()) != nullptr) {
                return false;
            }
        }
        return true;
    }() && "OptimizableSymmetryStorage::OptimizableSymmetryStorage: Cannot move from a SymmetryStorage that contains ReferenceSymmetryViews.");
}

OptimizableSymmetryStorage::~OptimizableSymmetryStorage() = default;

std::unique_ptr<ausaxs::symmetry::SymmetryStorage> OptimizableSymmetryStorage::clone() {
    auto copy = std::make_unique<OptimizableSymmetryStorage>(SymmetryStorage{});
    for (const auto& s : symmetries) { copy->symmetries.push_back(s->clone()); }
    return copy;
}
