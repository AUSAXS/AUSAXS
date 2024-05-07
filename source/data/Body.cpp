/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/state/UnboundSignaller.h>
#include <grid/Grid.h>
#include <constants/Constants.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <math/Vector3.h>

#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>

using namespace data;

Body::Body() {
    initialize();
}

Body::Body(const io::File& path) : uid(uid_counter++), file(path) {
    initialize();
}

Body::Body(const std::vector<record::Atom>& protein_atoms, const std::vector<record::Water>& hydration_atoms) : uid(uid_counter++), file(protein_atoms, hydration_atoms) {
    initialize();
}

Body::Body(const std::vector<record::Atom>& protein_atoms) : Body(protein_atoms, std::vector<record::Water>()) {}

Body::Body(const Body& body) : uid(body.uid), file(body.file) {
    initialize();
}

Body::Body(Body&& body) : uid(body.uid), file(std::move(body.file)) {
    initialize();
}

Body::~Body() = default;

void Body::initialize() {
    signal = std::make_shared<signaller::UnboundSignaller>();
}

void Body::save(const io::File& path) {file.write(path);}

void Body::center() {
    if (!centered) {
        translate(-get_cm());
        centered = true;
    }
}

Vector3<double> Body::get_cm() const {
    Vector3<double> cm;
    double M = 0; // total mass
    auto weighted_sum = [&cm, &M] (auto& atoms) {
        for (auto const& a : atoms) {
            double m = a.get_mass();
            M += m;
            cm += a.coords*m;
        }
    };
    weighted_sum(file.protein_atoms);
    weighted_sum(file.hydration_atoms);
    return cm/M;
}

double Body::get_volume_acids() const {
    double v = 0;
    int cur_seq = 0; // sequence number of current acid
    for (auto const& a : file.protein_atoms) {
        int a_seq = a.resSeq; // sequence number of current atom
        if (cur_seq != a_seq) { // check if we are still dealing with the same acid
            cur_seq = a_seq; // if not, update our current sequence number
            v += constants::volume::amino_acids.get(a.resName); // and add its volume to the running total
        }
    }
    return v;
}

void Body::translate(const Vector3<double>& v) {
    changed_external_state();

    std::for_each(file.protein_atoms.begin(), file.protein_atoms.end(), [&v] (record::Atom& atom) {atom.translate(v);});
    std::for_each(file.hydration_atoms.begin(), file.hydration_atoms.end(), [&v] (record::Water& atom) {atom.translate(v);});
}

void Body::rotate(const Matrix<double>& R) {
    for (auto& atom : file.protein_atoms) {
        atom.coords.rotate(R);
    }

    for (auto& atom : file.hydration_atoms) {
        atom.coords.rotate(R);
    }
}

void Body::rotate(double alpha, double beta, double gamma) {
    changed_external_state();
    Matrix R = matrix::rotation_matrix(alpha, beta, gamma);
    rotate(R);
}

void Body::rotate(const Vector3<double>& axis, double angle) {
    changed_external_state();
    Matrix R = matrix::rotation_matrix(axis, angle);
    rotate(R);
}

void Body::update_effective_charge(double charge) {
    changed_external_state();
    changed_internal_state();
    std::for_each(file.protein_atoms.begin(), file.protein_atoms.end(), [&charge] (record::Atom& a) {a.add_effective_charge(charge);});
    updated_charge = true;
}

double Body::get_total_atomic_charge() const {
    return std::accumulate(get_atoms().begin(), get_atoms().end(), 0.0, [] (double sum, const record::Atom& atom) {return sum + atom.Z();});
}

double Body::get_total_effective_charge() const {
    return std::accumulate(get_atoms().begin(), get_atoms().end(), 0.0, [](double sum, const record::Atom& a) { return sum + a.get_effective_charge(); });
}

double Body::get_molar_mass() const {
    return get_absolute_mass()*constants::Avogadro;
}

double Body::get_absolute_mass() const {
    double M = 0;
    std::for_each(file.protein_atoms.begin(), file.protein_atoms.end(), [&M] (const record::Atom& a) {M += a.get_mass();});
    std::for_each(file.hydration_atoms.begin(), file.hydration_atoms.end(), [&M] (const record::Water& a) {M += a.get_mass();});
    return M;
}

Body& Body::operator=(Body&& rhs) {
    file = std::move(rhs.file); 
    uid = rhs.uid;
    changed_internal_state();
    changed_external_state();
    return *this;
}

Body& Body::operator=(const Body& rhs) {
    file = rhs.file; 
    uid = rhs.uid;
    changed_internal_state();
    changed_external_state();
    return *this;
}

bool Body::operator==(const Body& rhs) const {
    return uid == rhs.uid;
}

bool Body::equals_content(const Body& rhs) const {
    return file.equals_content(rhs.file);
}

void Body::changed_external_state() const {signal->external_change();}

void Body::changed_internal_state() const {signal->internal_change();}

std::shared_ptr<signaller::Signaller> Body::get_signaller() const {
    return signal;
}

void Body::register_probe(std::shared_ptr<signaller::Signaller> signal) {
    this->signal = signal;
}

std::vector<record::Atom>& Body::get_atoms() {return file.protein_atoms;}

std::vector<record::Water>& Body::get_waters() {return file.hydration_atoms;}

const std::vector<record::Atom>& Body::get_atoms() const {return file.protein_atoms;}

const std::vector<record::Water>& Body::get_waters() const {return file.hydration_atoms;}

record::Atom& Body::get_atom(unsigned int index) {return file.protein_atoms[index];}

const record::Atom& Body::get_atom(unsigned int index) const {return file.protein_atoms[index];}

data::detail::AtomCollection& Body::get_file() {return file;}

int Body::get_id() const {return uid;}

std::size_t Body::size_atom() const {return file.protein_atoms.size();}