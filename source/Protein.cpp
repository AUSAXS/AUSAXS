#pragma once

// includes
#include <vector>

// my own includes
#include "Atom.cpp"
#include "pdbml_reader.h"

class Protein {
public:
    Protein(std::string filename) {
        Reader reader;
        if (filename.find(".xml") != std::string::npos) {
            reader = pdbml_reader(filename);
        }
        std::vector<Atom*> atoms = reader.read();
        separate(atoms);
    }

    /** Calculate the distances between each pair of atoms. 
     * @return A pair where the first entry is a vector of all internal distances between the protein atoms, while the second entry is all internal
     * distances between hydration atoms plus distances between hydration and protein atoms. 
     */
    std::pair<std::vector<double>, std::vector<double>> calc_distances() {
        // calculate the internal distances for the protein atoms
        int n = 0; // index counter
        std::vector<double> dp(protein.size()*(protein.size() - 1)); // n(n-1) total entries
        for (int i = 0; i < protein.size(); i++) {
            for (int j = i+1; j < protein.size(); j++) {
                dp[n] = protein[i]->distance(protein[j]);
                n++;
            }
        }
        std::cout << "n is: " << n << ", d is: " << dp.size() << std::endl;

        // calculate the distances for the hydrogen atoms
        n = 0; // index counter
        std::vector<double> dh(hydration.size()*(hydration.size() + protein.size() - 1)); // n(n-1) + nm = n(n + m - 1) total entries
        for (int i = 0; i < n; i++) {
            // loop over the hydration atoms
            for (int j = i+1; j < hydration.size(); j++) {
                dh[n] = hydration[i]->distance(hydration[j]);
                n++;
            }
            // loop over the protein atoms
            for (int j = 0; j < protein.size(); j++) {
                dh[n] = hydration[i]->distance(protein[j]);
                n++;
            }
        }
        std::cout << "n is: " << n << ", d is: " << dh.size() << std::endl;
        return std::make_pair(dp, dh);
    }
    std::vector<Atom*> protein; // atoms of the protein itself
    std::vector<Atom*> hydration; // hydration layer

private:

    /** Separate the structure into the protein and its hydration layer 
     * NOTE: consider removing return values
     * @return A pointer pair of (protein atoms, hydration atoms) to the private data of this class.
     */
    std::pair<std::vector<Atom*>*, std::vector<Atom*>*> separate(std::vector<Atom*> atoms) {
        hydration = std::vector<Atom*>(atoms.size());
        protein = std::vector<Atom*>(atoms.size());
        int i = 0, j = 0; // index counters for the hydration and protein vectors, respectively
        for (Atom* a : atoms) {
            if (a->get_comp() == "HOH") { // check if it is a hydration molecule
                hydration[i] = a;
                i++;
            } else {
                protein[j] = a;
                j++;
            }
        }
        hydration.resize(i);
        protein.resize(j);
        return std::make_pair(&protein, &hydration);
    };
};