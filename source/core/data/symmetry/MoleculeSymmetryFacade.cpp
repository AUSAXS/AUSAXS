// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <io/Writer.h>

#include <numeric>

using namespace ausaxs;
using namespace ausaxs::data;

std::size_t symmetry::detail::MoleculeSymmetryFacade::size_atom_total() const {
    return std::accumulate(
        molecule->get_bodies().begin(),
        molecule->get_bodies().end(),
        std::size_t{0},
        [] (std::size_t sum, const Body& body) {return sum + body.symmetry().size_atom_total();}
    );
}

std::size_t symmetry::detail::MoleculeSymmetryFacade::size_water_total() const {
    return std::accumulate(
        molecule->get_bodies().begin(),
        molecule->get_bodies().end(),
        std::size_t{0},
        [] (std::size_t sum, const Body& body) {return sum + body.symmetry().size_water_total();}
    );
}

data::detail::SimpleBody symmetry::detail::MoleculeSymmetryFacade::explicit_structure() const {
    std::vector<AtomFF> atoms;
    std::vector<Water> waters;
    atoms.reserve(size_atom_total());
    waters.reserve(size_water_total());
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