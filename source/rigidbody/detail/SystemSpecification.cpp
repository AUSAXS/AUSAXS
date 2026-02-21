// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/detail/SystemSpecification.h>

using namespace ausaxs::rigidbody::detail;

SystemSpecification::SystemSpecification() = default;

SystemSpecification::SystemSpecification(observer_ptr<const Rigidbody> rigidbody) 
    : absolute_parameters(rigidbody), initial_conformation(rigidbody->molecule.size_body())
{
    std::vector<data::Body> bodies = rigidbody->molecule.get_bodies();
    assert(absolute_parameters.parameters.size() == bodies.size() && "Configuration parameters size mismatch with molecule body size.");
    for (int i = 0; i < static_cast<int>(bodies.size()); ++i) {
        auto cm = bodies[i].get_cm();
        bodies[i].translate(-cm);
        initial_conformation[i] = std::move(bodies[i]);
        absolute_parameters.parameters[i].translation = cm;

        // initialize symmetry parameters from the body's symmetries
        for (int j = 0; j < static_cast<int>(initial_conformation[i].size_symmetry()); ++j) {
            absolute_parameters.parameters[i].symmetry_pars.emplace_back(initial_conformation[i].symmetry().get(j)->clone());
        }
    }
}

SystemSpecification::~SystemSpecification() = default;