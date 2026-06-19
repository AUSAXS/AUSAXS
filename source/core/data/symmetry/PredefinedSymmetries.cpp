// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/TetrahedralSymmetry.h>
#include <data/symmetry/OctahedralSymmetry.h>
#include <data/symmetry/IcosahedralSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>

#include <numbers>
#include <stdexcept>

std::unique_ptr<ausaxs::symmetry::ISymmetry> ausaxs::symmetry::get(type t) {
    switch (t) {
        case type::c2:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, std::numbers::pi},
                1
            );
        case type::c3:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/3},
                2
            );
        case type::c4:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/4},
                3
            );
        case type::c5:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/5},
                4
            );
        case type::c6:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/6},
                5
            );
        case type::p2:
            return std::make_unique<PointSymmetry>(
                Vector3<double>{0, 0, 0},
                Vector3<double>{0, 0, 0}
            );
        case type::t:
            return std::make_unique<TetrahedralSymmetry>();
        case type::o:
            return std::make_unique<OctahedralSymmetry>();
        case type::i:
            return std::make_unique<IcosahedralSymmetry>();
        default:
            throw std::runtime_error("Unknown symmetry type \"" + std::to_string(static_cast<int>(t)) + "\".");
    }
}

ausaxs::symmetry::type ausaxs::symmetry::get(std::string_view name) {
    if (name == "c2") {return type::c2;}
    if (name == "c3") {return type::c3;}
    if (name == "c4") {return type::c4;}
    if (name == "c5") {return type::c5;}
    if (name == "c6") {return type::c6;}
    if (name == "p2") {return type::p2;}
    if (name == "t")  {return type::t;}
    if (name == "o")  {return type::o;}
    if (name == "i")  {return type::i;}
    throw std::runtime_error("Unknown symmetry name \"" + std::string(name) + "\".");
}

std::unique_ptr<ausaxs::symmetry::ISymmetry> ausaxs::symmetry::create(std::string_view name) {
    // a hyphen denotes a nested composite: the first part is the inner symmetry, the remainder
    // (recursively parsed) the outer one
    if (auto pos = name.find('-'); pos != std::string_view::npos) {
        auto inner = create(name.substr(0, pos));
        auto outer = create(name.substr(pos + 1));
        return std::make_unique<CompositeSymmetry>(std::move(inner), std::move(outer));
    }
    return get(get(name));
}
