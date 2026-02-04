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

namespace {
    // Helper to compose rotation matrices and extract Euler angles (XYZ extrinsic convention)
    // The rotation matrix uses the convention: R = Rz(gamma) * Ry(beta) * Rx(alpha)
    ausaxs::Vector3<double> compose_rotation(const ausaxs::Matrix<double>& R_delta, const ausaxs::Vector3<double>& current_rotation) {
        auto R_current = ausaxs::matrix::rotation_matrix(current_rotation);
        auto R_new = R_delta * R_current;
        
        // Extract Euler angles (alpha, beta, gamma) from the rotation matrix
        double beta = std::asin(std::clamp(-R_new(2, 0), -1.0, 1.0));
        double cos_beta = std::cos(beta);
        
        double alpha, gamma;
        if (std::abs(cos_beta) > 1e-10) {
            alpha = std::atan2(R_new(2, 1), R_new(2, 2));
            gamma = std::atan2(R_new(1, 0), R_new(0, 0));
        } else {
            // Gimbal lock
            alpha = 0.0;
            gamma = std::atan2(-R_new(0, 1), R_new(1, 1));
        }
        
        return {alpha, beta, gamma};
    }
}

RigidTransform::RigidTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

RigidTransform::~RigidTransform() = default;

void RigidTransform::apply(parameter::BodyTransformParametersRelative&& par, constraints::DistanceConstraint& constraint) {
    auto group = get_connected(constraint);
    backup(group);

    auto grid = rigidbody->molecule.get_grid();

    // par is a delta transformation: par.rotation is the delta rotation, par.translation is the delta translation
    // The rotation happens around the group's pivot point.
    // For each body B with current absolute params (R_B, t_B):
    //   R_B_new = R_delta * R_B
    //   t_B_new = R_delta * (t_B - pivot) + pivot + t_delta
    auto R_delta = matrix::rotation_matrix(par.rotation);
    const auto& t_delta = par.translation;
    const auto& pivot = group.pivot;

    // Step 1: Compute new absolute parameters for all bodies in the group
    std::vector<std::pair<Vector3<double>, Vector3<double>>> new_params(group.bodies.size());
    
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        unsigned int ibody = group.indices[i];
        const auto& old_params = rigidbody->conformation->absolute_parameters.parameters[ibody];
        
        // R_B_new = R_delta * R_B
        // t_B_new = R_delta * (t_B - pivot) + pivot + t_delta
        auto new_rotation = compose_rotation(R_delta, old_params.rotation);
        auto new_translation = R_delta * (old_params.translation - pivot) + pivot + t_delta;
        
        new_params[i] = {new_rotation, new_translation};
    }
    
    // Step 2: Apply the computed absolute parameters to original_conformation
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        unsigned int ibody = group.indices[i];
        
        grid->remove(*group.bodies[i]);
        
        auto body = rigidbody->conformation->initial_conformation[ibody];
        auto& current_params = rigidbody->conformation->absolute_parameters.parameters[ibody];
        const auto& [rotation, translation] = new_params[i];
        body.rotate(matrix::rotation_matrix(rotation));
        body.translate(translation);
        
        current_params.rotation = rotation;
        current_params.translation = translation;

        auto& updated_body = rigidbody->molecule.get_body(ibody);
        current_params.symmetry_pars.clear();
        for (unsigned int i = 0; i < updated_body.size_symmetry(); ++i) {
            current_params.symmetry_pars.push_back(updated_body.symmetry().get(i));
        }
        
        rigidbody->molecule.get_body(ibody) = std::move(body);
    }

    // Handle symmetry for the first body
    if (!group.bodies.empty()) {
        symmetry(std::move(par.symmetry_pars), *group.bodies.front());
    }

    rigidbody->refresh_grid();

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