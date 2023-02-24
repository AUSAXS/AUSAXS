#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/RigidBody.h>

#include <vector>

using namespace rigidbody;

TransformStrategy::TransformGroup::TransformGroup(std::vector<Body*> bodies, std::vector<unsigned int> indices, std::shared_ptr<Constraint> target, Vector3<double> pivot) 
    : bodies(bodies), indices(indices), target(target), pivot(pivot) {}

TransformStrategy::TransformGroup TransformStrategy::get_connected(std::shared_ptr<Constraint> pivot) {
    // recursively explore a branch by stepping through its constraints, starting from the pivot and stopping if we reach ibody1 again
    unsigned int ibody1 = pivot->ibody1;
    std::function<std::vector<unsigned int>(unsigned int, std::vector<unsigned int>)> explore_branch = [&] (unsigned int ibody, std::vector<unsigned int> indices) {
        indices.push_back(ibody);
        if (ibody == ibody1) {
            return indices;
        }
        for (const auto& constraint : rigidbody->constraint_map[ibody]) {
            if (constraint->ibody1 == ibody) {
                explore_branch(constraint->ibody2, indices);
            } else {
                explore_branch(constraint->ibody1, indices);
            }
        }
        return indices;
    };

    // explore all branches
    auto path1 = explore_branch(pivot->ibody1, std::vector<unsigned int>());
    auto path2 = explore_branch(pivot->ibody2, std::vector<unsigned int>());

    // if the paths are the same length, we just return the pivot as the only body in the group
    if (path1.size() == path2.size()) {
        return TransformGroup({&rigidbody->bodies[ibody1]}, {ibody1}, pivot, pivot->get_atom1().coords);
    }

    // create a vector of pointers to the bodies in the paths
    std::vector<Body*> bodies1, bodies2;
    for (const auto& ibody : path1) {
        bodies1.push_back(&rigidbody->bodies[ibody]);
    }
    for (const auto& ibody : path2) {
        bodies2.push_back(&rigidbody->bodies[ibody]);
    }

    // check if the system is overconstrained
    if (0.5*rigidbody->body_size() < path1.size() && 0.5*rigidbody->body_size() < path2.size()) {
        throw except::size_error("TransformStrategy::get_connected: The system is overconstrained. Use a different TransformStrategy.");
    }

    // if the paths are different lengths, we return the shorter path as the group
    if (path1.size() < path2.size()) {
        return TransformGroup(bodies1, path1, pivot, pivot->get_atom1().coords);
    } else {
        return TransformGroup(bodies2, path2, pivot, pivot->get_atom2().coords);
    }
}

void TransformStrategy::rotate(const Matrix<double>& M, TransformGroup& group) {
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (Body* body) {body->translate(-group.pivot);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&M] (Body* body) {body->rotate(M);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (Body* body) {body->translate(group.pivot);});
}

void TransformStrategy::translate(const Vector3<double>& t, TransformGroup& group) {
    std::for_each(group.bodies.begin(), group.bodies.end(), [&t] (Body* body) {body->translate(t);});
}

void TransformStrategy::undo() {
    for (auto& body : bodybackup) {
        rigidbody->bodies[body.index] = std::move(body.body);
    }
    bodybackup.clear();
}

void TransformStrategy::backup(TransformGroup& group) {
    bodybackup.clear();
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        bodybackup.emplace_back(*group.bodies[i], group.indices[i]);
    }
}
