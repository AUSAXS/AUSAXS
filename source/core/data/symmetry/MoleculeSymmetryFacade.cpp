#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/Writer.h>

using namespace ausaxs;
using namespace ausaxs::data;

std::vector<AtomFF> symmetry::detail::MoleculeSymmetryFacade::explicit_structure() const {
    std::vector<AtomFF> atoms = molecule->get_atoms();
    int N = std::accumulate(
        molecule->get_bodies().begin(), 
        molecule->get_bodies().end(), 
        0, 
        [] (int sum, const Body& body) {return sum + body.size_atom()*(body.size_symmetry_total()+1);}
    );
    atoms.reserve(N);
    for (const auto& body : molecule->get_bodies()) {
        auto body_atoms = body.symmetry().explicit_structure().get_atoms();
        atoms.insert(atoms.end(), body_atoms.begin(), body_atoms.end());
    }
    return atoms;
}

void symmetry::detail::MoleculeSymmetryFacade::save(const io::File& path) const {
    auto atoms = explicit_structure();
    io::Writer::write(io::pdb::PDBStructure(Body(std::move(atoms))), path);
}