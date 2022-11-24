#include <data/BodySplitter.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <utility/Exceptions.h>

#include <algorithm>

Protein BodySplitter::split(const std::string input, std::vector<int> splits) {
    Body body(input);
    std::vector<Atom>& atoms = body.atoms();

    // we define a boolean vector with one entry for each residue sequence id
    std::vector<bool> split(atoms.back().resSeq, false);

    // we then mark the ids where we want to split as true
    std::for_each(splits.begin(), splits.end(), [&split] (const int id) {split[id] = true;});

    std::vector<Body> bodies(splits.size()+1);
    int index_body = 0; // current index in the bodies vector

    // the two iterators marks the indices in the atoms vector where we want to split next time
    std::vector<Atom>::const_iterator begin = atoms.begin(); // start at the beginning
    std::vector<Atom>::const_iterator end; // no defined end yet
    for (unsigned int i = 0; i < atoms.size(); i++) {
        int resSeq = std::max(atoms[i].resSeq, 0); // in some files resSeq starts negative

        // we can now in constant time look in our split vector to see if we should split at this atom
        if (split[resSeq]) {
            end = atoms.begin() + i;        // define the end index
            std::vector<Atom> a(begin, end);     // create a new vector of atoms based on the start and end iterators
            bodies[index_body++] = Body(a); // create a body from this vector
            begin = end;                    // change the start index for the next split
            split[resSeq] = false;          // mark it as false so we won't split again on the next atom
        }
    }

    // add the final body
    begin = end;
    end = atoms.end();
    std::vector<Atom> a(begin, end);
    bodies[index_body] = Body(a);

    return Protein(bodies);
}

std::vector<Constraint> BodySplitter::sequential_constraints(const Protein& protein) {
    const std::vector<Body>& bodies = protein.bodies;

    std::vector<Constraint> constraints(bodies.size()-1);
    for (unsigned int i = 0; i < protein.bodies.size()-1; i++) {
        const Body &body1 = bodies[i], &body2 = bodies[i+1]; 

        int res1 = body1.atoms().back().resSeq;
        int res2 = body2.atoms()[0].resSeq;

        unsigned int index1 = -1, index2 = -1;
        for (unsigned int j = body1.atoms().size()-1; j > 0; j--) {
            const Atom& atom = body1.atoms(j);
            if (atom.resSeq != res1) [[unlikely]] { // sanity check
                throw except::unexpected("BodySplitter::sequential_constrain: Could not find C-alpha atom.");
            }
            if (atom.name == "CA") {
                index1 = j;
                break;
            }
        }
        for (unsigned int j = 0; j < body2.atoms().size(); j++) {
            const Atom& atom = body2.atoms(j);
            if (atom.resSeq != res2) [[unlikely]] { // sanity check
                throw except::unexpected("BodySplitter::sequential_constrain: Could not find C-alpha atom.");
            }
            if (atom.name == "CA") {
                index2 = j;
                break;
            }
        }

        // sanity check
        if (res1 == -1 || res2 == -1) {
            throw except::unexpected("BodySplitter::sequential_constrain: Could not find C-alpha atom.");
        }
        constraints[i] = Constraint(&body1.atoms(index1), &body2.atoms(index2), &body1, &body2);
    }
    return constraints;
}