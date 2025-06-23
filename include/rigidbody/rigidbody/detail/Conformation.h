#pragma once

#include <rigidbody/detail/Configuration.h>
#include <data/Body.h>

namespace ausaxs::rigidbody::detail {
    /**
     * @brief The current conformation of a rigid-body molecule.
     *        Storing the original conformation is required to express the current conformation in absolute terms. 
     */
    struct Conformation {
        Conformation() = default;
        Conformation(const std::vector<data::Body>& bodies) : configuration(), original_conformation(bodies) {}

        Configuration configuration;
        std::vector<data::Body> original_conformation;
    };
}