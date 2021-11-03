#pragma once

// includes
#include <vector>
#include <map>
#include "boost/format.hpp"

// ROOT
#include <TVector3.h>

// my own includes
#include "Atom.cpp"
#include "pdbml_reader.cpp"

using boost::format;
using std::vector, std::string, std::cout, std::endl;
using namespace ROOT;

class Protein {
public:
    Protein(string filename) {
        Reader* reader;
        if (filename.find(".xml") != string::npos) {
            reader = new pdbml_reader(filename);
        }
        vector<Atom*> atoms = reader->read();
        separate(atoms);
    }

    /** Calculate the distances between each pair of atoms. 
     * @return A pair where the first entry is a vector of all internal distances between the protein atoms, while the second entry is all internal
     * distances between hydration atoms plus distances between hydration and protein atoms. 
     */
    std::pair<vector<double>, vector<double>> calc_distances() {
        // calculate the internal distances for the protein atoms
        int n = 0; // index counter
        vector<double> dp(protein_atoms.size()*(protein_atoms.size() - 1)/2); // n(n-1)/2 total entries
        for (int i = 0; i < protein_atoms.size(); i++) {
            for (int j = i+1; j < protein_atoms.size(); j++) {
                dp[n] = protein_atoms[i]->distance(protein_atoms[j]);
                n++;
            }
        }

        // calculate the distances for the hydrogen atoms
        n = 0; // index counter
        vector<double> dh(hydration_atoms.size()*(hydration_atoms.size() + 2*protein_atoms.size() - 1)/2); // n(n-1)/2 + nm = n(n + 2m - 1)/2 total entries
        for (int i = 0; i < hydration_atoms.size(); i++) {
            // loop over the hydration atoms
            for (int j = i+1; j < hydration_atoms.size(); j++) {
                dh[n] = hydration_atoms[i]->distance(hydration_atoms[j]);
                n++;
            }
            // loop over the protein atoms
            for (int j = 0; j < protein_atoms.size(); j++) {
                dh[n] = hydration_atoms[i]->distance(protein_atoms[j]);
                n++;
            }
        }
        return make_pair(dp, dh);
    }

    /** Use an algorithm to generate a new hydration layer for this protein. Note that the previous one will be deleted.
     * 
     */
    void generate_new_hydration() {
        // delete the old hydration layer
        for (Atom* a : hydration_atoms) {
            delete a;
        }

        // move protein to centre of mass
        TVector3 cm = get_cm();
        translate(-cm);

        // generate the 3D grid
        double width = 10; // what width to use? 10 is too large, but with smaller values our grid becomes incredibly large
        auto[corner, bins] = gridify(width); // corner is the lower corner of our grid, and bins the number of bins in each dimension
        cout << format("bins: (%1%, %2%, %3%)") % bins[0] % bins[1] % bins[2] << endl;
        vector<vector<vector<bool>>> occupied = find_protein_locations(corner, bins, width);

        int n = 0;
        for (int i = 0; i < bins[0]; i++) {
            for (int j = 0; j < bins[1]; j++) {
                for (int k = 0; k < bins[2]; k++) {
                    if (occupied[i][j][k]) {
                        n++;
                    }
                }
            }
        }
        cout << format("Occupied slots: %1%, number of protein atoms: %2%") % n % protein_atoms.size() << endl;
    }

private:
    vector<Atom*> protein_atoms; // atoms of the protein itself
    vector<Atom*> hydration_atoms; // hydration layer

    /** Marks the approximate protein locations (and *not* hydration locations!) in the input grid. Helper function for generate_new_hydration. 
     * @return A boolean array which is "true" if the grid point is occupied, and "false" otherwise. 
     */
    vector<vector<vector<bool>>> find_protein_locations(const TVector3 corner, const vector<int> bins, const double width) {
        vector<vector<vector<bool>>> occupied(bins[0], vector<vector<bool>>(bins[1], vector<bool>(bins[2], false))); // THIS ARRAY IS HUUGE!!!
        for (Atom* a : protein_atoms) {
            cout << "Hello there atom no " << a->get_serial() << endl;
            occupied.at(int((a->get_x() - corner.X())/width)-1).at(int((a->get_y() - corner.Y())/width)-1).at(int((a->get_z() - corner.Z())/width)-1) = true;
            // occupied[int((a->get_x() - corner.X())/width)][int((a->get_y() - corner.Y())/width)][int((a->get_y() - corner.Y())/width)] = true;
        }
        return occupied;
    }

    /** Generate a 3D grid containing all atoms. Helper function for generate_new_hydration. 
     * @param width the bin width
     * @return A pair (corner point, bins) where the first is a bottom corner of the grid, while the second is the number of bins in each dimension.
     */
    std::pair<TVector3, vector<int>> gridify(const double width) {
        // determine the size of our grid
        TVector3 high = get_cm();
        TVector3 low = high;
        auto update = [&low, &high] (vector<Atom*>* atoms) {
            for (Atom* a : *atoms) {
                // update minimum vector
                if (a->get_x() < low.X()) low.SetX(a->get_x());
                if (a->get_y() < low.Y()) low.SetY(a->get_y());
                if (a->get_z() < low.Z()) low.SetZ(a->get_z());

                // update maximum vector
                if (a->get_x() > high.X()) high.SetX(a->get_x());
                if (a->get_y() > high.Y()) high.SetY(a->get_y());
                if (a->get_z() > high.Z()) high.SetZ(a->get_z());
            }
        };
        update(&protein_atoms);
        update(&hydration_atoms);

        // calculate the number of bins in each dimension and initialize the grid and occupancy vectors
        vector<int> bins = {int((high.X() - low.X())/width), int((high.Y() - low.Y())/width), int((high.Z() - low.Z())/width)};
        return std::make_pair(low, bins);
    }

    /** Move the entire protein by a vector
     * @param v the translation vector
     */
    void translate(const TVector3 v) {
        auto move = [&v] (vector<Atom*>* atoms) {
            for (Atom* a : *atoms) {
                a->translate(v);
            }
        };
        move(&protein_atoms);
        move(&hydration_atoms);
    }

    /** Calculate the center-mass coordinates for the protein.
     * @return The center-mass (x, y, z) coordinates. 
     */
    TVector3 get_cm() {
        TVector3 cm;
        auto weighted_sum = [&cm] (vector<Atom*>* atoms) {
            for (Atom* a : *atoms) {
                double w = a->get_atomic_weight();
                double x = a->get_x()/w;
                double y = a->get_y()/w;
                double z = a->get_z()/w;
                cm += TVector3(x, y, z);
            }
        };
        weighted_sum(&protein_atoms);
        weighted_sum(&hydration_atoms);
        return cm;
    }

    /** Separate the structure into the protein and its hydration layer 
     * NOTE: consider removing return values
     * @return A pointer pair of (protein atoms, hydration atoms) to the private data of this class.
     */
    std::pair<vector<Atom*>*, vector<Atom*>*> separate(const vector<Atom*> atoms) {
        hydration_atoms = vector<Atom*>(atoms.size());
        protein_atoms = vector<Atom*>(atoms.size());
        int i = 0, j = 0; // index counters for the hydration and protein vectors, respectively
        for (Atom* a : atoms) {
            if (a->get_comp() == "HOH") { // check if it is a hydration molecule
                hydration_atoms[i] = a;
                i++;
            } else {
                protein_atoms[j] = a;
                j++;
            }
        }
        hydration_atoms.resize(i);
        protein_atoms.resize(j);
        cout << format("Found %1% protein atoms, and %2% hydration atoms.") % protein_atoms.size() % hydration_atoms.size() << endl;
        return make_pair(&protein_atoms, &hydration_atoms);
    };
};