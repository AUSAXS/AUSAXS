#pragma once

#include <rigidbody/detail/Configuration.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>

namespace ausaxs::rigidbody::detail {
    /**
     * @brief The current conformation of a rigid-body molecule.
     *        Storing the original conformation is required to express the current conformation in absolute terms. 
     */
    struct Conformation {
        Conformation();
        Conformation(observer_ptr<const Rigidbody> rigidbody);
        ~Conformation();

        Configuration configuration;
        std::vector<data::Body> original_conformation;
    };
}