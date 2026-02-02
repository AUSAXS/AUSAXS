// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/parameters/BodyTransformParameters.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/Rigidbody.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/Body.h>

#include <unordered_set>
#include <functional>
#include <numeric>

using namespace ausaxs::rigidbody::transform;

RigidTransform::RigidTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

RigidTransform::~RigidTransform() = default;

void RigidTransform::apply(parameter::BodyTransformParameters&& par, constraints::DistanceConstraint& constraint) {
    auto group = get_connected(constraint);
    backup(group);

    auto grid = rigidbody->molecule.get_grid();

    // Calculate transformation delta for the primary body
    auto& old_primary_params = rigidbody->conformation->configuration.parameters[constraint.ibody1];
    auto delta_rotation = par.rotation - old_primary_params.rotation;
    auto delta_rotation_matrix = matrix::rotation_matrix(delta_rotation);

    // Get the current center of mass of the primary body to use as rotation pivot
    auto primary_center = rigidbody->molecule.get_body(constraint.ibody1).get_cm();

    // For each body in the rigid group, apply absolute transformation from original_conformation
    // The primary body gets the new parameters, secondary bodies maintain rigid relationship
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        unsigned int ibody = group.indices[i];

        // Remove old body from grid
        grid->remove(*group.bodies[i]);

        // Get fresh body from original_conformation
        auto body = rigidbody->conformation->original_conformation[ibody];

        // Apply absolute transformation using parameters
        if (ibody == constraint.ibody1) {
            // Primary body: use the new parameters
            body.rotate(matrix::rotation_matrix(par.rotation));
            body.translate(par.translation);
            symmetry(std::move(par.symmetry_pars), body);

            // Update configuration with new absolute parameters
            rigidbody->conformation->configuration.parameters[ibody] = std::move(par);
        } else {
            // Secondary bodies: same rotation delta, but translation accounts for orbital motion
            auto& current_params = rigidbody->conformation->configuration.parameters[ibody];
            
            // Apply same rotation delta to maintain relative orientation
            auto new_rotation = current_params.rotation + delta_rotation;
            
            // Calculate new translation: body orbits around primary's center during rotation
            auto current_center = rigidbody->molecule.get_body(ibody).get_cm();
            auto offset_from_primary = current_center - primary_center;
            auto rotated_offset = delta_rotation_matrix * offset_from_primary;
            auto orbital_displacement = rotated_offset - offset_from_primary;
            auto new_translation = current_params.translation + (par.translation - old_primary_params.translation) + orbital_displacement;

            body.rotate(matrix::rotation_matrix(new_rotation));
            body.translate(new_translation);

            // Update configuration with new absolute parameters
            rigidbody->conformation->configuration.parameters[ibody].rotation = new_rotation;
            rigidbody->conformation->configuration.parameters[ibody].translation = new_translation;
        }

        // Update molecule with transformed body
        rigidbody->molecule.get_body(ibody) = std::move(body);
    }

    // Ensure grid has space for new conformations
    rigidbody->refresh_grid();

    // Add all bodies back to grid
    for (auto& body : group.bodies) {
        grid->add(*body);
    }
}

TransformGroup RigidTransform::get_connected(const constraints::DistanceConstraint& pivot) {
    // explore the graph of bodies connected to 'ibody'
    std::function<void(unsigned int, std::unordered_set<unsigned int>&)> explore_branch = [&] (unsigned int ibody, std::unordered_set<unsigned int>& indices) {
        // if we've already explored this branch, return
        if (indices.contains(ibody)) {
            return;
        }
        // otherwise, add this body to the list of explored bodies
        indices.insert(ibody);

        // explore all bodies connected to this body
        for (const auto& constraint : rigidbody->constraints->distance_constraints_map[ibody]) {
            if (constraint.get().ibody1 == ibody) {
                explore_branch(constraint.get().ibody2, indices);
            } else {
                explore_branch(constraint.get().ibody1, indices);
            }
        }
        return;
    };

    // explore all branches
    std::unordered_set<unsigned int> _path1({pivot.ibody2});
    std::unordered_set<unsigned int> _path2({pivot.ibody1});
    explore_branch(pivot.ibody1, _path1);
    explore_branch(pivot.ibody2, _path2);
    _path1.erase(pivot.ibody2);
    _path2.erase(pivot.ibody1);
    std::vector<unsigned int> path1(_path1.begin(), _path1.end());
    std::vector<unsigned int> path2(_path2.begin(), _path2.end());

    // if the paths are the same length, we just return the pivot as the only body in the group
    if (path1.size() == path2.size() && path1 == path2) {
        return TransformGroup({&rigidbody->molecule.get_body(pivot.ibody1)}, {pivot.ibody1}, pivot, pivot.get_atom1().coordinates());
    }

    // create a vector of pointers to the bodies in the paths
    std::vector<observer_ptr<data::Body>> bodies1, bodies2;
    for (const auto& ibody : path1) {
        bodies1.push_back(&rigidbody->molecule.get_body(ibody));
    }
    for (const auto& ibody : path2) {
        bodies2.push_back(&rigidbody->molecule.get_body(ibody));
    }

    // check if the system is overconstrained
    if (0.5*rigidbody->molecule.size_body() < path1.size() && 0.5*rigidbody->molecule.size_body() < path2.size()) {
        throw except::size_error("TransformStrategy::get_connected: The system is overconstrained. Use a different TransformStrategy.");
    }

    unsigned int N1 = std::accumulate(path1.begin(), path1.end(), 0, [&] (unsigned int sum, unsigned int ibody) {
        return sum + rigidbody->molecule.get_body(ibody).size_atom();
    });
    unsigned int N2 = std::accumulate(path2.begin(), path2.end(), 0, [&] (unsigned int sum, unsigned int ibody) {
        return sum + rigidbody->molecule.get_body(ibody).size_atom();
    });

    // return the path with the least atoms, since that will be the cheapest to transform
    if (N1 < N2) {
        return TransformGroup(std::move(bodies1), std::move(path1), pivot, pivot.get_atom2().coordinates());
    } else {
        return TransformGroup(std::move(bodies2), std::move(path2), pivot, pivot.get_atom1().coordinates());
    }
}