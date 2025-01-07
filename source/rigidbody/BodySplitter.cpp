/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/BodySplitter.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <io/pdb/PDBStructure.h>
#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <io/Reader.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::data;
using namespace ausaxs::io::pdb;

Molecule BodySplitter::split(const io::File& input, std::vector<int> splits) {
    io::pdb::PDBStructure data = io::Reader::read(input);
    std::vector<PDBAtom>& atoms = data.atoms;

    // we define a boolean vector with one entry for each residue sequence id
    int max_id = 0;
    for (const auto& a : atoms) {max_id = std::max(max_id, a.resSeq);}
    std::vector<bool> split_at(max_id, false);

    // we then mark the ids where we want to split as true
    std::for_each(splits.begin(), splits.end(), [&split_at] (unsigned int id) {
        if (id > split_at.size()) {throw except::parse_error("BodySplitter::split: Split index (" + std::to_string(id) + ") larger than highest residue sequence id (" + std::to_string(split_at.size()) + ").");}
        split_at[id] = true;
    });

    std::vector<Body> bodies(splits.size()+1);
    int index_body = 0; // current index in the bodies vector

    // the two iterators marks the indices in the atoms vector where we want to split next time
    auto begin = atoms.begin();
    for (unsigned int i = 0; i < atoms.size(); i++) {
        int resSeq = std::max(atoms[i].resSeq, 0); // in some files resSeq starts negative

        if (split_at[resSeq]) {
            std::vector<PDBAtom> a(begin, atoms.begin()+i);
            auto reduced = PDBStructure(a, {}).reduced_representation();
            bodies[index_body++] = Body(reduced.atoms, reduced.waters);
            split_at[resSeq] = false; // mark it as false so we won't split again on the next atom
            begin = atoms.begin() + i;
        }
    }
    auto reduced = PDBStructure(std::vector<PDBAtom>(begin, atoms.end()), {}).reduced_representation();
    bodies[index_body] = Body(reduced.atoms, reduced.waters);
    return Molecule(bodies);
}

data::Molecule BodySplitter::split(const io::File& input) {
    io::pdb::PDBStructure data = io::Reader::read(input);
    std::vector<PDBAtom>& atoms = data.atoms;

    std::vector<Body> bodies;
    auto begin = atoms.begin();
    char current_id = atoms[0].chainID;
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (atoms[i].chainID != current_id) {
            std::vector<PDBAtom> a(begin, atoms.begin() + i);
            auto reduced = PDBStructure(a, {}).reduced_representation();
            bodies.push_back(Body(reduced.atoms, reduced.waters));
            begin = atoms.begin() + i;
            current_id = atoms[i].chainID;
        }
    }
    auto reduced = PDBStructure(std::vector<PDBAtom>(begin, atoms.end()), {}).reduced_representation();
    bodies.push_back(Body(reduced.atoms, reduced.waters));
    return Molecule(bodies);
}