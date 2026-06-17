#pragma once

#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <data/symmetry/ReferenceSymmetry.h>

#include <cassert>

using namespace ausaxs::symmetry;

OptimizableSymmetryStorage::OptimizableSymmetryStorage(SymmetryStorage&& other) : SymmetryStorage(std::move(other)) {
    assert([] () {
        for (const auto& s : symmetries) {
            if (dynamic_cast<const ReferenceSymmetryView*>(s.get()) != nullptr) {
                return false;
            }
        }
        return true;
    } && "OptimizableSymmetryStorage::OptimizableSymmetryStorage: Cannot move from a SymmetryStorage that contains ReferenceSymmetryViews.");
}
OptimizableSymmetryStorage::~OptimizableSymmetryStorage() = default;

inline std::unique_ptr<ausaxs::symmetry::SymmetryStorage> ausaxs::symmetry::OptimizableSymmetryStorage::clone() {
    auto copy = std::make_unique<ausaxs::symmetry::OptimizableSymmetryStorage>(ausaxs::symmetry::SymmetryStorage{});
    for (const auto& s : symmetries) { copy->symmetries.push_back(s->clone()); }
    copy->optimize_translate = optimize_translate;
    copy->optimize_rot_axis = optimize_rot_axis;
    return copy;
}

inline void ausaxs::symmetry::OptimizableSymmetryStorage::add(symmetry::type symmetry) {
    symmetries.emplace_back(symmetry::get(symmetry));
    switch (symmetry) {
        case symmetry::type::c2:
        case symmetry::type::c3:
        case symmetry::type::c4:
        case symmetry::type::c5:
        case symmetry::type::c6:
        case symmetry::type::p2:
        case symmetry::type::t:
        case symmetry::type::o:
        case symmetry::type::i:
            optimize_rot_axis = true;
            optimize_translate = true;
            break;

        default:
            throw std::runtime_error("OptimizableSymmetryStorage::add: Symmetry \"" + std::to_string(static_cast<int>(symmetry)) + "\" not found. Did you forget to add it?");
    }
}