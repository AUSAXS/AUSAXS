#include <rigidbody/parameters/Parameters.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/Exceptions.h>
#include <math/Vector3.h>
#include <Symbols.h>

using namespace rigidbody::parameter;

Parameters::Parameters(const data::Molecule* protein) : params(protein->body_size()) {
    const std::vector<data::Body>& bodies = protein->get_bodies();
    for (unsigned int i = 0; i < params.size(); i++) {
        id_to_index[bodies[i].get_id()] = i;
    }
}

void Parameters::update(unsigned int uid, const Parameter& param) {
    #ifdef DEBUG
        if (id_to_index.find(uid) == id_to_index.end()) {throw except::out_of_bounds("Parameters::update: Body uid \"" + std::to_string(uid) + "\" not found");}
    #endif

    params[id_to_index[uid]] = param;
}

void Parameters::update(unsigned int uid, Vector3<double> dx, double drx, double dry, double drz) {
    update(uid, Parameter(dx, drx, dry, drz));
}

const Parameter Parameters::get(unsigned int uid) {
    return params[id_to_index[uid]];
}

std::string Parameter::to_string() const {
    return "translation: " + dr.to_string() + ", angles: (" + std::to_string(alpha) + ", " + std::to_string(beta) + ", " + std::to_string(gamma) + ")"; 
}