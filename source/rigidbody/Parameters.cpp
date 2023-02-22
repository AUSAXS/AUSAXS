#include <rigidbody/parameters/Parameters.h>

#include <vector>

using namespace rigidbody;

Parameters::Parameters(const Protein& protein) : params(protein.bodies.size()) {
    const std::vector<Body>& bodies = protein.bodies;
    for (unsigned int i = 0; i < params.size(); i++) {
        id_to_index[bodies[i].uid] = i;
    }
}

void Parameters::update(unsigned int uid, const Parameter& param) {
    params[id_to_index[uid]] = param;
}

void Parameters::update(unsigned int uid, Vector3<double> dx, double drx, double dry, double drz) {
    update(uid, Parameter(dx, drx, dry, drz));
}

const Parameter Parameters::get(unsigned int uid) {
    return params[id_to_index[uid]];
}