#pragma once

#include <rigidbody/detail/Configuration.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>

namespace ausaxs::rigidbody::detail {
    /**
     * @brief The current conformation of a rigid-body molecule.
     * 
     * This structure holds all original bodies and their current configuration.
     *
     * For convenience, the contained bodies are moved to the origin, with the configuration
     * initialized with the center of mass of each body required to restore the original conformation.
     */
    struct Conformation {
        Conformation();
        Conformation(observer_ptr<const Rigidbody> rigidbody);
        ~Conformation();

        Configuration configuration;
        std::vector<data::Body> original_conformation;
    };
}