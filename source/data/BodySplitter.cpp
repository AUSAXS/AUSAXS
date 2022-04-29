#pragma once

#include <data/BodySplitter.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <Exceptions.h>

#include <algorithm>

using std::vector, std::string;

Protein BodySplitter::split(const string input, vector<int> splits) {
    Body body(input);
    vector<Atom>& atoms = body.protein_atoms;

    // we define a boolean vector with one entry for each residue sequence id
    vector<bool> split(atoms.back().resSeq, false);

    // we then mark the ids where we want to split as true
    std::for_each(splits.begin(), splits.end(), [&split] (const int id) {split[id] = true;});

    vector<Body> bodies(splits.size()+1);
    int index_body = 0; // current index in the bodies vector

    // the two iterators marks the indices in the atoms vector where we want to split next time
    vector<Atom>::const_iterator begin = atoms.begin(); // start at the beginning
    vector<Atom>::const_iterator end; // no defined end yet
    for (unsigned int i = 0; i < atoms.size(); i++) {
        int resSeq = std::max(atoms[i].resSeq, 0); // in some files resSeq starts negative

        // we can now in constant time look in our split vector to see if we should split at this atom
        if (split[resSeq]) {
            end = atoms.begin() + i;        // define the end index
            vector<Atom> a(begin, end);     // create a new vector of atoms based on the start and end iterators
            bodies[index_body++] = Body(a); // create a body from this vector
            begin = end;                    // change the start index for the next split
            split[resSeq] = false;          // mark it as false so we won't split again on the next atom
        }
    }

    // add the final body
    begin = end;
    end = atoms.end();
    vector<Atom> a(begin, end);
    bodies[index_body] = Body(a);

    return Protein(bodies);
}

vector<Constraint> BodySplitter::sequential_constraints(const Protein& protein) {
    const vector<Body>& bodies = protein.bodies;

    vector<Constraint> constraints(bodies.size()-1);
    for (unsigned int i = 0; i < protein.bodies.size()-1; i++) {
        const Body &body1 = bodies[i], &body2 = bodies[i+1]; 

        int res1 = body1.protein_atoms.back().resSeq;
        int res2 = body2.protein_atoms[0].resSeq;

        unsigned int index1 = -1, index2 = -1;
        for (unsigned int j = body1.protein_atoms.size()-1; j > 0; j--) {
            const Atom& atom = body1.protein_atoms[j];
            if (__builtin_expect(atom.resSeq != res1, false)) { // sanity check
                throw except::unexpected("Error in BodySplitter::sequential_constrain: Could not find C-alpha atom.");
            }
            if (atom.name == "CA") {
                index1 = j;
                break;
            }
        }
        for (unsigned int j = 0; j < body2.protein_atoms.size(); j++) {
            const Atom& atom = body2.protein_atoms[j];
            if (__builtin_expect(atom.resSeq != res2, false)) { // sanity check
                throw except::unexpected("Error in BodySplitter::sequential_constrain: Could not find C-alpha atom.");
            }
            if (atom.name == "CA") {
                index2 = j;
                break;
            }
        }

        // sanity check
        if (res1 == -1 || res2 == -1) {
            throw except::unexpected("Error in BodySplitter::sequential_constrain: Could not find C-alpha atom.");
        }
        constraints[i] = Constraint(&body1.protein_atoms[index1], &body2.protein_atoms[index2], &body1, &body2);
    }
    return constraints;
}