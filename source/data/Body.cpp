#include "data/Body.h"

Body::Body(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms) : file(FileConstructor::construct()), protein_atoms(file->protein_atoms), hydration_atoms(file->hydration_atoms){}
Body::Body(const string& input) : file(FileConstructor::construct(input)), protein_atoms(file->protein_atoms), hydration_atoms(file->hydration_atoms) {
    file->read();
}

void Body::translate(const Vector3& v) {
    auto move = [&v] (auto& atoms) {
        for (auto& a : atoms) {
            a.translate(v);
        }
    };
    move(protein_atoms);
    move(hydration_atoms);
}
void Body::rotate(const double& alpha, const double& beta, const double& gamma) {
    double sina = 
}