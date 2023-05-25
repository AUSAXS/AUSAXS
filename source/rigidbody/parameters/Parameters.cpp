#include <rigidbody/parameters/Parameters.h>
#include <data/Protein.h>
#include <data/Body.h>

#include <vector>

using namespace rigidbody;

Parameter::Parameter() : dx(0, 0, 0), alpha(0), beta(0), gamma(0) {}

Parameter::Parameter(const Vector3<double>& dx, double alpha, double beta, double gamma) : dx(dx), alpha(alpha), beta(beta), gamma(gamma) {}

Parameters::Parameters(const Protein* protein) : params(protein->body_size()) {
    const std::vector<Body>& bodies = protein->get_bodies();
    for (unsigned int i = 0; i < params.size(); i++) {
        id_to_index[bodies[i].get_id()] = i;
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

std::string Parameter::to_string() const {
    return "translation: " + dx.to_string() + ", angles: (" + std::to_string(alpha) + ", " + std::to_string(beta) + ", " + std::to_string(gamma) + ")"; 
}