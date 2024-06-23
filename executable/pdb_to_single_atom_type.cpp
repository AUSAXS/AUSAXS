#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>

int main(int argc, char const *argv[]) {
    io::ExistingFile pdb = argv[1];
    data::Molecule data(pdb);

    for (auto& body : data.get_bodies()) {
        for (auto& atom : body.get_atoms()) {
            atom.set_element(constants::atom_t::C);
            atom.set_residue_name("UNK");
            atom.set_group_name("C");
        }
    }

    data.save(pdb.directory().path() + "/carbonified.pdb");
    return 0;
}