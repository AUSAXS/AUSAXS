#pragma once

#include <data/DataFwd.h>
#include <data/symmetry/Symmetry.h>
#include <utility/observer_ptr.h>

namespace ausaxs::symmetry {
    enum type {
        p2,
        p3,
        p4
    };

    /**
     * @brief Apply a given symmetry to a body.
     */
    void apply(observer_ptr<data::Body> body, type t);
}