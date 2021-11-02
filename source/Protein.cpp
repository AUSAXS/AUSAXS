#pragma once

// includes
#include <vector>

// my own includes
#include <Atom.cpp>

class Protein {
public:
    std::vector<Atom*> atoms; // all atoms
    std::vector<Atom*> protein; // atoms of the protein itself
    std::vector<Atom*> hydration; // hydration layer

    Protein(std::vector<Atom*> atoms) {
        this->atoms = atoms;
    }

private:
    /** Separate the structure into the hydration layer and the protein itself
     * @return a pair of (protein atoms, hydration atoms)
     */
    std::pair<std::vector<Atom*>, std::vector<Atom*>> separate() {};

};