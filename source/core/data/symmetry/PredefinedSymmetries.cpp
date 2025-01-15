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