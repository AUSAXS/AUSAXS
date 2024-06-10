/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/BodySplitter.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <rigidbody/constraints/Constraint.h>

#include <algorithm>

using namespace rigidbody;
using namespace data;
using namespace data::record;

Molecule BodySplitter::split(const io::File& input, std::vector<int> splits) {
    Body body(input);
    std::vector<Atom>& atoms = body.get_atoms();

    // we define a boolean vector with one entry for each residue sequence id
    int max_id = 0;
    for (const auto& a : atoms) {max_id = std::max(max_id, a.resSeq);}
    std::vector<bool> split_at(max_id, false);

    // we then mark the ids where we want to split as true
    std::for_each(splits.begin(), splits.end(), [&split_at] (unsigned int id) 
        {
            if (id > split_at.size()) {throw except::parse_error("BodySplitter::split: Split index (" + std::to_string(id) + ") larger than highest residue sequence id (" + std::to_string(split_at.size()) + ").");}
            split_at[id] = true;
        }
    );

    std::vector<Body> bodies(splits.size()+1);
    int index_body = 0; // current index in the bodies vector

    // the two iterators marks the indices in the atoms vector where we want to split next time
    std::vector<Atom>::const_iterator begin = atoms.begin(); // start at the beginning
    std::vector<Atom>::const_iterator end; // no defined end yet
    for (unsigned int i = 0; i < atoms.size(); i++) {
        int resSeq = std::max(atoms[i].resSeq, 0); // in some files resSeq starts negative

        // we can now in constant time look in our split vector to see if we should split at this atom
        if (split_at[resSeq]) {
            end = atoms.begin() + i;        // define the end index
            std::vector<Atom> a(begin, end);// create a new vector of atoms based on the start and end iterators
            bodies[index_body++] = Body(a); // create a body from this vector
            begin = end;                    // change the start index for the next split
            split_at[resSeq] = false;       // mark it as false so we won't split again on the next atom
        }
    }

    // add the final body
    begin = end;
    end = atoms.end();
    std::vector<Atom> a(begin, end);
    bodies[index_body] = Body(a);

    return Molecule(bodies);
}