#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/RigidBody.h>

#include <vector>

using namespace rigidbody;

// TransformStrategy::TransformGroup TransformStrategy::get_connected(const Constraint& pivot) const {
//     const std::vector<Constraint>& constraints = rigidbody->constraints;    // easy access to the set of all constraints
//     std::unordered_map<unsigned int, unsigned int> group;                   // a map of body uids to group ids
//     std::unordered_map<unsigned int, std::list<const Body*>> group_members; // a map of group ids to a list of its members

//     // each body starts in its own group
//     std::for_each(rigidbody->bodies.begin(), rigidbody->bodies.end(), [&group, &group_members] (const Body& body) {
//         group[body.uid] = body.uid;
//         group_members[body.uid] = {&body};
//     });

//     // iterate through all constraints
//     for (unsigned int i = 0; i < constraints.size(); i++) {
//         const Constraint& constraint = constraints[i];

//         // if the constraint is the pivot, we skip it
//         if (constraint == pivot) {
//             continue;
//         }

//         // get the current group id of each body from the constraint
//         unsigned int id1 = group[constraint.get_body1().uid];
//         unsigned int id2 = group[constraint.get_body2().uid];

//         // if they are not already in the same group, we merge their groups
//         if (group[id1] != group[id2]) {
//             // references to the two group lists
//             std::list<const Body*>& members1 = group_members[id1];
//             std::list<const Body*>& members2 = group_members[id2];

//             // change all members of group2 to group1
//             for (const auto& body : members2) {
//                 group[body->uid] = id1;
//                 members1.push_back(body);
//             }
//             group_members.erase(id2); // erase the old member list
//         }
//     }

//     // get the group id of the pivot, and return a vector of all bodies from the same group
//     unsigned int id = group[pivot.get_body1().uid];
//     std::vector<Body*> connected;
//     for(auto& body : rigidbody->bodies) {
//         if (group[body.uid] == id) {
//             connected.push_back(&body);
//         }
//     }
//     return TransformGroup{.bodies=connected, .pivot=&pivot};
// }

TransformStrategy::TransformGroup TransformStrategy::get_connected(const Constraint& pivot) {
    // check if the constraint map is up to date
    if (constraint_map.size() != rigidbody->bodies.size()) {
        generate_constraint_map();
    }

    // recursively explore a branch by stepping through its constraints, starting from the pivot and stopping if we reach ibody1 again
    unsigned int ibody1 = pivot.ibody1;
    std::function<std::vector<Body*>(unsigned int, std::vector<Body*>)> explore_branch = [&] (unsigned int ibody, std::vector<Body*> branch) {
        branch.push_back(&rigidbody->bodies[ibody]);
        if (ibody == ibody1) {
            return branch;
        }
        for (const auto& constraint : constraint_map[ibody]) {
            if (constraint->ibody1 == ibody) {
                explore_branch(constraint->ibody2, branch);
            } else {
                explore_branch(constraint->ibody1, branch);
            }
        }
        return branch;
    };

    // explore all branches
    auto path1 = explore_branch(pivot.ibody1, std::vector<Body*>());
    auto path2 = explore_branch(pivot.ibody2, std::vector<Body*>());

    // if the paths are the same length, we just return the pivot as the only body in the group
    if (path1.size() == path2.size()) {
        return TransformGroup{.bodies={&rigidbody->bodies[ibody1]}, .pivot=&pivot};
    }

    if (0.5*rigidbody->body_size() < path1.size() && 0.5*rigidbody->body_size() < path2.size()) {
        throw except::size_error("TransformStrategy::get_connected: The system is overconstrained. Use a different TransformStrategy.");
    }

    // if the paths are different lengths, we return the shorter path as the group
    if (path1.size() < path2.size()) {
        return TransformGroup{.bodies=path1, .pivot=&pivot};
    } else {
        return TransformGroup{.bodies=path2, .pivot=&pivot};
    }
}

void TransformStrategy::generate_constraint_map() {
    for (const auto& body : rigidbody->bodies) {
        constraint_map[body.uid] = std::vector<std::shared_ptr<rigidbody::Constraint>>();
    }

    for (const auto& constraint : rigidbody->get_constraints()) {
        constraint_map[constraint->get_body1().uid].push_back(constraint);
        constraint_map[constraint->get_body2().uid].push_back(constraint);
    }
}