#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/Writer.h>

using namespace ausaxs;
using namespace ausaxs::data;

data::detail::SimpleBody symmetry::detail::MoleculeSymmetryFacade::explicit_structure() const {
    std::vector<AtomFF> atoms;
    std::vector<Water> waters;
    int Na = std::accumulate(
        molecule->get_bodies().begin(), 
        molecule->get_bodies().end(), 
        0, 
        [] (int sum, const Body& body) {return sum + body.symmetry().size_atom_total();}
    );
    int Nw = std::accumulate(
        molecule->get_bodies().begin(), 
        molecule->get_bodies().end(), 
        0, 
        [] (int sum, const Body& body) {return sum + body.symmetry().size_water_total();}
    );
    atoms.reserve(Na);
    waters.reserve(Nw);
    for (const auto& body : molecule->get_bodies()) {
        auto body_atoms = body.symmetry().explicit_structure();
        atoms.insert(atoms.end(), body_atoms.atoms.begin(), body_atoms.atoms.end());
        waters.insert(waters.end(), body_atoms.waters.begin(), body_atoms.waters.end());
    }
    return {atoms, waters};
}

bool symmetry::detail::MoleculeSymmetryFacade::has_symmetries() const {
    return std::accumulate(
        molecule->get_bodies().begin(), 
        molecule->get_bodies().end(), 
        false, 
        [] (bool sum, const Body& body) {return sum || body.size_symmetry();}
    );
}

void symmetry::detail::MoleculeSymmetryFacade::save(const io::File& path) const {
    auto body = explicit_structure();
    io::Writer::write(io::pdb::PDBStructure(Body(std::move(body.atoms), std::move(body.waters))), path);
}