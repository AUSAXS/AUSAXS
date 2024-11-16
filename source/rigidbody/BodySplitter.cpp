/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/BodySplitter.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <rigidbody/constraints/Constraint.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::data;
using namespace ausaxs::data::record;

Molecule BodySplitter::split(const io::File& input, std::vector<int> splits) {
    Body body(input);
    std::vector<Atom>& atoms = body.get_atoms();

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
            std::vector<Atom> a(begin, atoms.begin()+i);
            bodies[index_body++] = Body(a);
            split_at[resSeq] = false; // mark it as false so we won't split again on the next atom
            begin = atoms.begin() + i;
        }
    }
    bodies[index_body] = Body(std::vector<Atom>(begin, atoms.end()));
    return Molecule(bodies);
}

data::Molecule BodySplitter::split(const io::File& input) {
    Body body(input);
    std::vector<Atom>& atoms = body.get_atoms();

    std::vector<Body> bodies;
    auto begin = atoms.begin();
    char current_id = atoms[0].chainID;
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (atoms[i].chainID != current_id) {
            std::vector<Atom> a(begin, atoms.begin() + i);
            bodies.push_back(Body(a));
            begin = atoms.begin() + i;
            current_id = atoms[i].chainID;
        }
    }
    bodies.push_back(Body(std::vector<Atom>(begin, atoms.end())));
    return Molecule(bodies);
}