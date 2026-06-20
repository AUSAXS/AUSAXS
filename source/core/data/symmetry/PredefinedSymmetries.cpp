// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/TetrahedralSymmetry.h>
#include <data/symmetry/OctahedralSymmetry.h>
#include <data/symmetry/IcosahedralSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <utility/StringUtils.h>

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
        case type::c7:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/7},
                6
            );
        case type::c8:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/8},
                7
            );
        case type::c9:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/9},
                8
            );
        case type::c10:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/10},
                9
            );
        case type::c11:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/11},
                10
            );
        case type::c12:
            return std::make_unique<CyclicSymmetry>(
                CyclicSymmetry::_Relation{{0, 0, 0}},
                CyclicSymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/12},
                11
            );
        case type::p2:
            return std::make_unique<PointSymmetry>(
                Vector3<double>{0, 0, 0},
                Vector3<double>{0, 0, 0}
            );
        case type::tetrahedral:
            return std::make_unique<TetrahedralSymmetry>();
        case type::octahedral:
            return std::make_unique<OctahedralSymmetry>();
        case type::icosahedral:
            return std::make_unique<IcosahedralSymmetry>();
        default:
            throw std::runtime_error("Unknown symmetry type \"" + std::to_string(static_cast<int>(t)) + "\".");
    }
}

ausaxs::symmetry::type ausaxs::symmetry::get(std::string_view name) {
    if (name.empty()) {throw std::runtime_error("symmetry::get: Symmetry name cannot be empty.");}
    if (name == "c2") {return type::c2;}
    if (name == "c3") {return type::c3;}
    if (name == "c4") {return type::c4;}
    if (name == "c5") {return type::c5;}
    if (name == "c6") {return type::c6;}
    if (name == "c7") {return type::c7;}
    if (name == "c8") {return type::c8;}
    if (name == "c9") {return type::c9;}
    if (name == "c10") {return type::c10;}
    if (name == "c11") {return type::c11;}
    if (name == "c12") {return type::c12;}
    if (name == "p2") {return type::p2;}
    if (name == "t" || name == "tetrahedral") {return type::tetrahedral;}
    if (name == "o" || name == "octahedral")  {return type::octahedral;}
    if (name == "i" || name == "icosahedral") {return type::icosahedral;}
    throw std::runtime_error("symmetry::get: Unknown symmetry name \"" + std::string(name) + "\".");
}

std::unique_ptr<ausaxs::symmetry::ISymmetry> ausaxs::symmetry::create(std::string_view name) {
    std::string lc = utility::to_lowercase(name);
    // a hyphen denotes a nested composite: the first part is the inner symmetry, the remainder (recursively parsed) the outer one
    if (auto pos = lc.find('-'); pos != std::string_view::npos) {
        auto inner = create(lc.substr(0, pos));
        auto outer = create(lc.substr(pos + 1));
        return std::make_unique<CompositeSymmetry>(std::move(inner), std::move(outer));
    }
    return get(get(lc));
}
