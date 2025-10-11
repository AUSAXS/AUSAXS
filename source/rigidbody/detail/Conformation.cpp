// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/detail/Conformation.h>

using namespace ausaxs::rigidbody::detail;

Conformation::Conformation() = default;

Conformation::Conformation(observer_ptr<const Rigidbody> rigidbody) 
    : configuration(rigidbody), original_conformation(rigidbody->molecule.size_body())
{
    std::vector<data::Body> bodies = rigidbody->molecule.get_bodies();
    assert(configuration.parameters.size() == bodies.size() && "Configuration parameters size mismatch with molecule body size.");
    for (unsigned int i = 0; i < bodies.size(); ++i) {
        auto cm = bodies[i].get_cm();
        bodies[i].translate(-cm);
        original_conformation[i] = std::move(bodies[i]);
        configuration.parameters[i].translation = cm;
    }
}

Conformation::~Conformation() = default;