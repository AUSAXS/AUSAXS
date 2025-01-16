#pragma once

#include <data/DataFwd.h>
#include <data/symmetry/Symmetry.h>
#include <utility/observer_ptr.h>

namespace ausaxs::symmetry {
    enum class type {
        p2,
        p3,
        p4
    };

    /**
     * @brief Apply a given symmetry to a body.
     */
    symmetry::Symmetry get(type t);

    /**
     * @brief Get a predefined symmetry by name.
     */
    type get(std::string_view name);
}