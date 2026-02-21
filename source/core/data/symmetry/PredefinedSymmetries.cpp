// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/Symmetry.h>
#include <data/symmetry/PointSymmetry.h>

#include <numbers>
#include <stdexcept>

std::unique_ptr<ausaxs::symmetry::ISymmetry> ausaxs::symmetry::get(type t) {
    switch (t) {
        case type::c2:
            return std::make_unique<Symmetry>(
                ISymmetry::_Relation{{0, 0, 0}},
                ISymmetry::_Repeat{{0, 0, 1}, std::numbers::pi},
                1
            );
        case type::c3:
            return std::make_unique<Symmetry>(
                ISymmetry::_Relation{{0, 0, 0}},
                ISymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/3},
                2
            );
        case type::c4:
            return std::make_unique<Symmetry>(
                ISymmetry::_Relation{{0, 0, 0}},
                ISymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/4},
                3
            );
        case type::c5:
            return std::make_unique<Symmetry>(
                ISymmetry::_Relation{{0, 0, 0}},
                ISymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/5},
                4
            );
        case type::c6:
            return std::make_unique<Symmetry>(
                ISymmetry::_Relation{{0, 0, 0}},
                ISymmetry::_Repeat{{0, 0, 1}, 2*std::numbers::pi/6},
                5
            );
        case type::p2:
            // Decoupled dimer: zero initial offset and orientation along z-axis.
            // The optimiser will freely perturb translation (d) and orientation (axis/angle).
            return std::make_unique<PointSymmetry>(
                Vector3<double>{0, 0, 0},
                Vector3<double>{0, 0, 1},
                std::numbers::pi
            );
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
    throw std::runtime_error("Unknown symmetry name \"" + std::string(name) + "\".");
}
