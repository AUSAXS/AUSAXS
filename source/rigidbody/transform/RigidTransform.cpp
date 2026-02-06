// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>

#include <unordered_set>
#include <functional>
#include <numeric>

using namespace ausaxs::rigidbody::transform;

RigidTransform::RigidTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

RigidTransform::~RigidTransform() = default;

void RigidTransform::apply(parameter::BodyTransformParametersRelative&& par, constraints::DistanceConstraint& constraint) {
    auto group = get_connected(constraint);
    backup(group);
    auto grid = rigidbody->molecule.get_grid();

    // Step 1: Get fresh bodies and compute new absolute parameters
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        unsigned int ibody = group.indices[i];
        auto& body = *group.bodies[i];
        grid->remove(body);
        body = rigidbody->conformation->initial_conformation[ibody];
    }

    // Step 2: Compute new absolute parameters for all bodies
    std::vector<parameter::BodyTransformParametersAbsolute> new_params(group.bodies.size());
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        unsigned int ibody = group.indices[i];
        new_params[i] = rigidbody->conformation->absolute_parameters.parameters[ibody];
        
        if (par.rotation.has_value()) {
            new_params[i].rotation += par.rotation.value();
        }
        if (par.translation.has_value()) {
            auto R_delta = matrix::rotation_matrix(par.rotation.value_or(Vector3<double>{0, 0, 0}));
            new_params[i].translation = R_delta * (new_params[i].translation - group.pivot) + group.pivot + par.translation.value();
        }
        if (i == 0 && par.symmetry_pars.has_value()) {
            new_params[i].symmetry_pars = add_symmetries(new_params[i].symmetry_pars, par.symmetry_pars.value());
        }
    }
    
    // Step 3: Apply transformations to all bodies
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        auto& body = *group.bodies[i];
        const auto& params = new_params[i];
        
        body.rotate(matrix::rotation_matrix(params.rotation));
        body.translate(params.translation);
        
        if (i == 0 && par.symmetry_pars.has_value()) {
            apply_symmetry(params.symmetry_pars, body);
        }
        
        grid->add(body);
    }

    rigidbody->refresh_grid();
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