#include <data/symmetry/PredefinedSymmetries.h>
#include <data/Body.h>

ausaxs::symmetry::Symmetry ausaxs::symmetry::get(type t) {
    switch (t) {
        case type::p2:
            return {
                {
                    {0, 0, 0},
                    {0, 0, 0}
                },
                {
                    {0, 0, 0},
                    {0, 0, std::numbers::pi/2}
                }
            };
            break;

        case type::p3:
            return {
                {
                    {0, 0, 0},
                    {0, 0, 0}
                },
                {
                    {0, 0, 0},
                    {0, 0, std::numbers::pi/3}
                }
            };
            break;

        case type::p4:
            return {
                {
                    {0, 0, 0},
                    {0, 0, 0}
                },
                {
                    {0, 0, 0},
                    {0, 0, std::numbers::pi/4}
                }
            };
            break;

        default: 
            throw std::runtime_error("Unknown symmetry type \"" + std::to_string(static_cast<int>(t)) + "\".");
    }
}

ausaxs::symmetry::type ausaxs::symmetry::get(std::string_view name) {
    if (name == "p2") {return type::p2;}
    if (name == "p3") {return type::p3;}
    if (name == "p4") {return type::p4;}
    throw std::runtime_error("Unknown symmetry name \"" + std::string(name) + "\".");
}