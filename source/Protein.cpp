#pragma once

// includes
#include <vector>
#include <map>
#include "boost/format.hpp"

// ROOT
#include <TVector3.h>

// my own includes
#include "Atom.cpp"
#include "io/PDBML_reader.cpp"
#include "io/PDB_reader.cpp"
#include "io/PDBML_writer.cpp"
#include "io/PDB_writer.cpp"

using boost::format;
using std::vector, std::string, std::cout, std::endl;
using namespace ROOT;

class Protein {
public:
    /** Creates a new protein from the input .pdb or .xml file. 
     * @param path path to the input file. 
     */
    Protein(string path) {
        // determine which kind of input file we're looking at
        Reader* reader;
        if (path.find(".xml") != string::npos) { // .xml file
            reader = new PDBML_reader(path);
        } else if (path.find(".pdb") != string::npos) { // .pdb file
            reader = new PDB_reader(path);
        } else { // anything else - we cannot handle this
            print_err((format("ERROR: Invalid file extension of input file %1%.") % path).str());
            exit(1);
        }
        vector<Atom*> atoms = reader->read();
        reader->close();

        separate(atoms);
    }

    /** Writes this protein to disk.
     * @param path path to the destination. 
     */
    void save(string path) {
        // determine the format of the output file
        Writer* writer;
        if (path.find(".xml") != string::npos) { // .xml file
            writer = new PDBML_writer(path);
        } else if (path.find(".pdb") != string::npos) { // .pdb file
            writer = new PDB_writer(path);
        } else { // anything else - we cannot handle this
            print_err((format("ERROR: Unknown writing format for path \"%1%\".") % path).str());
            exit(1);
        }
        writer->write(&protein_atoms, &hydration_atoms);
        writer->close();
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
        hydration_atoms = vector<Atom*>();

        // move protein to centre of mass
        TVector3 cm = get_cm();
        cout << format("Center-of-mass is (x, y, z) = (%1%, %2%, %3%)") % cm[0] % cm[1] % cm[2] << endl;
        translate(-cm);

        // generate the 3D grid
        double width = 10; // what width to use? 10 is too large, but with smaller values our grid becomes incredibly large
        auto[corner, bins] = generate_grid(width); // corner is the lower corner of our grid, and bins the number of bins in each dimension
        cout << format("bins: (%1%, %2%, %3%)") % bins[0] % bins[1] % bins[2] << endl;
        // vector<vector<vector<bool>>> occupied = find_protein_locations(corner, bins, width);
        

        // Jans approach: generate grid from -250 to 250 in all dimensions, with grid size 1Å. 
        // Move protein to center so its more likely to be covered by the grid. 
        // Determine which grid points are filled by our atoms. We can do this in a smart way by converting each coordinate to an index instead of searching through the entire grid. 
        // Expand each grid point into a proper size for its atom. Jan uses sqrt(6) in each direction. 
        // Calculate the volume by iterating through the grid and counting how many entries are filled. Can probably be done smarter by keeping track of it along the way.
        //      The volume may be scaled. Jan multiplied by 27 for some reason. 
        // Determine the surface atoms of the protein. Jan did this by checking 3Å in each direction. If no other atom occupies the spot, it is a candidate for a hydration atom. 
        // Place the hydration atoms. Only keep every third. 
    }

    vector<Atom*>* get_protein_atoms() {
        return &protein_atoms;
    }

    vector<Atom*>* get_hydration_atoms() {
        return &hydration_atoms;
    }

    /** Generate a 3D grid containing all atoms.
     * @param width the bin width
     * @return A pair (corner point, bins) where the first is a bottom corner of the grid, while the second is the number of bins in each dimension.
     */
    std::pair<TVector3, vector<int>> generate_grid(const double width) {
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

    /** Calculate the center-mass coordinates for the protein.
     * @return The center-mass (x, y, z) coordinates. 
     */
    TVector3 get_cm() {
        TVector3 cm;
        double M = 0; // total mass
        auto weighted_sum = [&cm, &M] (vector<Atom*>* atoms) {
            for (Atom* a : *atoms) {
                double m = a->get_atomic_weight();
                M += m;
                double x = a->get_x()*m;
                double y = a->get_y()*m;
                double z = a->get_z()*m;
                cm += TVector3(x, y, z);
            }
            cm[0] = cm[0]/M;
            cm[1] = cm[1]/M;
            cm[2] = cm[2]/M;
        };
        weighted_sum(&protein_atoms);
        weighted_sum(&hydration_atoms);
        return cm;
    }

    /** Generate a boolean grid array representation of the protein. Occupied slots are marked by "true". 
     * @param locs the index locations of each atom.
     * @return A 3D boolean array representation of the protein. 
     */
    // vector<vector<vector<bool>>> boolean_representation(std::map<int, vector<int>> locs) {
        // double width = 0.1;
        // double radius = 2.5; // radius of each atom in Ångström
        // auto[corner, bins] = generate_grid(width);
        // std::map<int, vector<int>> locs = find_protein_locations(corner, bins, width);

        // cout << format("Number of bins (x, y, z): (%1%, %2%, %3%)") % bins[0] % bins[1] % bins[2] << endl;

        // // create the boolean grid array
        // vector<vector<vector<bool>>> grid(bins[0], vector<vector<bool>>(bins[1], vector<bool>(bins[2], false)));

        // return grid;
    // }

private:
    vector<Atom*> protein_atoms; // atoms of the protein itself
    vector<Atom*> hydration_atoms; // hydration layer

    void place_hydration_atoms() {
        // determine the surface
        // vector<int> hydration_locs();
        // auto check_dim = [&] (int dim) {
        //     if (i-radius_bins > 0) {
        //         if (!grid[i-radius_bins][j][k]) {
        //             hydration_locs.push_back({i - radius_bins, j, k});
        //         }
        //     }
        // };
        // for (int i = 0; i < bins[0]; i++) {
        //     for (int j = 0; j < bins[1]; j++) {
        //         for (int k = 0; k < bins[2]; k++) {
        //             if (!grid[i][j][k]) continue; // empty

        //             // now we just check all possible locations

        //         }
        //     }
        // }        
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