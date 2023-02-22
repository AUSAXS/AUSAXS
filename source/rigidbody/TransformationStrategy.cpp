#include <rigidbody/transform/TransformationStrategy.h>
#include <rigidbody/RigidBody.h>

#include <vector>

using namespace rigidbody;

std::vector<Body*> TransformationStrategy::get_connected(const Constraint& pivot) const {
    const std::vector<Constraint>& constraints = protein->constraints;  // easy access to the set of all constraints
    std::unordered_map<size_t, size_t> group;                           // a map of body uids to group ids
    std::unordered_map<size_t, std::list<const Body*>> group_members;   // a map of group ids to a list of its members

    // each body starts in its own group
    std::for_each(protein->protein.bodies.begin(), protein->protein.bodies.end(), [&group, &group_members] (const Body& body) {
        group[body.uid] = body.uid;
        group_members[body.uid] = {&body};
    });

    // iterate through all constraints
    for (unsigned int i = 0; i < constraints.size(); i++) {
        const Constraint& constraint = constraints[i];

        // if the constraint is the pivot, we skip it
        if (constraint == pivot) {
            continue;
        }

        // get the current group id of each body from the constraint
        unsigned int id1 = group.at(constraint.get_body1().uid);
        unsigned int id2 = group.at(constraint.get_body2().uid);

        // if they are not already in the same group, we merge their groups
        if (group.at(id1) != group.at(id2)) {
            // references to the two group lists
            std::list<const Body*>& members1 = group_members.at(id1);
            std::list<const Body*>& members2 = group_members.at(id2);

            // change all members of group2 to group1
            for (const auto& body : members2) {
                group[body->uid] = id1;
                members1.push_back(body);
            }
            group_members.erase(id2); // erase the old member list
        }
    }

    // get the id of the pivot, and return a vector of all bodies from the same group
    unsigned int id = group.at(pivot.get_body1().uid);
    std::vector<Body*> connected;
    for(auto& body : protein->protein.bodies) {
        if (group.at(body.uid) == id) {
            connected.push_back(&body);
        }
    }
    return connected;
}