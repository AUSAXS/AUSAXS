#include "data/Protein.h"
#include "data/Body.h"

Protein::Protein(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms) : file(std::make_unique<File>()) {}
Protein::Protein(const string& input) : file(std::make_unique<File>(input)) {}

void Protein::translate(const Vector3& v) {
    for (auto& body : bodies) {
        body.translate(v);
    }
}