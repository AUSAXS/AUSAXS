// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PredefinedSymmetries.h>
#include <data/Body.h>

ausaxs::symmetry::Symmetry ausaxs::symmetry::get(type t) {
    switch (t) {
        case type::c2:
            return {
                {{0, 0, 0}},
                {{0, 0, 1}, std::numbers::pi},
                1
            };

        case type::c3:
            return {
                {{0, 0, 0}},
                {{0, 0, 1}, 2*std::numbers::pi/3},
                2
            };

        case type::c4:
            return {
                {{0, 0, 0}},
                {{0, 0, 1}, 2*std::numbers::pi/4},
                3
            };
            break;

        case type::c5:
            return {
                {{0, 0, 0}},
                {{0, 0, 1}, 2*std::numbers::pi/5},
                4
            };

        case type::c6:
            return {
                {{0, 0, 0}},
                {{0, 0, 1}, 2*std::numbers::pi/6},
                5
            };

        case type::p2:
            return {
                {{0, 0, 0}},
                {{0, 0, 1}, std::numbers::pi},
                1
            };

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