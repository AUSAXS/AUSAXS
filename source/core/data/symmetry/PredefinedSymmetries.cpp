#include <data/symmetry/PredefinedSymmetries.h>
#include <data/Body.h>

void ausaxs::symmetry::apply(observer_ptr<data::Body> body, type t) {
    switch (t) {
        case p2:
            body->symmetry().add({
                {0, 0, 0},
                {0, 0, 0},
                {0, 0, std::numbers::pi/2}
            });
            break;

        case p3:
            body->symmetry().add({
                {0, 0, 0},
                {0, 0, 0},
                {0, 0, std::numbers::pi/3}
            });
            break;

        case p4:
            body->symmetry().add({
                {0, 0, 0},
                {0, 0, 0},
                {0, 0, std::numbers::pi/4}
            });
            break;

        default: 
            throw std::runtime_error("Unknown symmetry type \"" + std::to_string(t) + "\".");
    }
}

ausaxs::symmetry::type ausaxs::symmetry::get(std::string_view name) {
    if (name == "p2") {return p2;}
    if (name == "p3") {return p3;}
    if (name == "p4") {return p4;}
    throw std::runtime_error("Unknown symmetry name \"" + std::string(name) + "\".");
}