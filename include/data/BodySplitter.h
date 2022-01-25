#pragma once

#include <algorithm>

#include "data/Body.h"
#include "data/Protein.h"
#include "data/Atom.h"

using std::vector, std::string;

struct BodySplitter {
    /**
     * @brief Constructor. 
     * 
     * @param input The path to the input data file. 
     */
    static Protein split(const string input, vector<int> splits) {
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
                end = atoms.begin() + i;                              // define the end index
                vector<Atom> a(begin, end);                           // create a new vector of atoms based on the start and end iterators
                bodies[index_body++] = Body(a);                       // create a body from this vector
                begin = end;                                          // change the start index for the next split
                split[resSeq] = false;                                // mark it as false so we won't split again on the next atom
            }
        }

        // add the final body
        begin = end;
        end = atoms.end();
        vector<Atom> a(begin, end);
        bodies[index_body] = Body(a);

        return Protein(bodies);
    }
};