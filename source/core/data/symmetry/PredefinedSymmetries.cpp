// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/TetrahedralSymmetry.h>
#include <data/symmetry/OctahedralSymmetry.h>
#include <data/symmetry/IcosahedralSymmetry.h>
#include <data/symmetry/DihedralSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <utility/StringUtils.h>

#include <algorithm>
#include <cctype>
#include <numbers>
#include <optional>
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
        case type::d2:  return std::make_unique<DihedralSymmetry>(2);
        case type::d3:  return std::make_unique<DihedralSymmetry>(3);
        case type::d4:  return std::make_unique<DihedralSymmetry>(4);
        case type::d5:  return std::make_unique<DihedralSymmetry>(5);
        case type::d6:  return std::make_unique<DihedralSymmetry>(6);
        case type::d7:  return std::make_unique<DihedralSymmetry>(7);
        case type::d8:  return std::make_unique<DihedralSymmetry>(8);
        case type::d9:  return std::make_unique<DihedralSymmetry>(9);
        case type::d10: return std::make_unique<DihedralSymmetry>(10);
        case type::d11: return std::make_unique<DihedralSymmetry>(11);
        case type::d12: return std::make_unique<DihedralSymmetry>(12);
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
    if (name == "d2")  {return type::d2;}
    if (name == "d3")  {return type::d3;}
    if (name == "d4")  {return type::d4;}
    if (name == "d5")  {return type::d5;}
    if (name == "d6")  {return type::d6;}
    if (name == "d7")  {return type::d7;}
    if (name == "d8")  {return type::d8;}
    if (name == "d9")  {return type::d9;}
    if (name == "d10") {return type::d10;}
    if (name == "d11") {return type::d11;}
    if (name == "d12") {return type::d12;}
    throw std::runtime_error("symmetry::get: Unknown symmetry name \"" + std::string(name) + "\".");
}

namespace {
    // parse a bare cyclic token "c<k>" into its order k; std::nullopt if it is not one
    std::optional<int> parse_cyclic(const std::string& s) {
        if (s.size() < 2 || s[0] != 'c') {return std::nullopt;}
        if (!std::all_of(s.begin() + 1, s.end(), [](unsigned char c) {return std::isdigit(c);})) {return std::nullopt;}
        return std::stoi(s.substr(1));
    }

    // "c2" combined with a single "c<n>" (either order) shares one centre with perpendicular axes,
    // which is the dihedral group D_n; returns n when the two tokens match that pattern
    std::optional<int> dihedral_order(const std::string& left, const std::string& right) {
        auto a = parse_cyclic(left), b = parse_cyclic(right);
        if (!a || !b || std::min(*a, *b) != 2) {return std::nullopt;}
        return std::max(*a, *b);
    }
}

std::unique_ptr<ausaxs::symmetry::ISymmetry> ausaxs::symmetry::create(std::string_view name) {
    std::string lc = utility::to_lowercase(name);
    // a hyphen denotes a nested composite: the first part is the inner symmetry, the remainder (recursively parsed) the outer one
    if (auto pos = lc.find('-'); pos != std::string_view::npos) {
        std::string left = lc.substr(0, pos), right = lc.substr(pos + 1);
        // ... except a bare "c2-cN" / "cN-c2" pair is the dihedral group D_N, built explicitly so
        // it forms a genuine 2N-copy point group and its pair schedule exploits the full symmetry
        if (auto n = dihedral_order(left, right)) {return get(get("d" + std::to_string(*n)));}
        return std::make_unique<CompositeSymmetry>(create(left), create(right));
    }
    return get(get(lc));
}
