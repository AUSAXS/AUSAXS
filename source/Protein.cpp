#pragma once

// includes
#include <vector>
#include <map>
#include "boost/format.hpp"
#include <utility>

// ROOT
#include <TVector3.h>

// my own includes
#include "data/Atom.cpp"
#include "Grid.cpp"
#include "data/PDB_file.cpp"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::unique_ptr;
using namespace ROOT;

class Protein {
public:
    /** Creates a new protein from the input .pdb or .xml file. 
     * @param path path to the input file. 
     */
    Protein(string path) {
        // determine which kind of input file we're looking at
        if (path.find(".xml") != string::npos) { // .xml file
            print_err("Error in Protein::Protein: .xml input files are not supported.");
        } else if (path.find(".pdb") != string::npos) { // .pdb file
            file = std::make_shared<PDB_file>(path);
        } else { // anything else - we cannot handle this
            print_err((format("Error in Protein::Protein: Invalid file extension of input file %1%.") % path).str());
            exit(1);
        }
        
        std::tie(protein_atoms, hydration_atoms) = file->get_atoms();
    }

    /** Writes this protein to disk.
     * @param path path to the destination. 
     */
    void save(string path) {
        file->update(protein_atoms, hydration_atoms); // update the File backing this Protein with our new atoms
        file->write(path); // write to disk
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

    /** 
     * @brief Use an algorithm to generate a new hydration layer for this protein. Note that the previous one will be deleted.
     */
    void generate_new_hydration() {
        // delete the old hydration layer
        hydration_atoms = vector<shared_ptr<Atom>>();

        // move protein to center of mass
        TVector3 cm = get_cm();
        cout << format("Center-of-mass is (x, y, z) = (%1%, %2%, %3%)") % cm[0] % cm[1] % cm[2] << endl;
        translate(-cm);

        // generate the 3D grid
        Grid grid({-250, -250, -250}, 1, 501);
        grid.add(&protein_atoms);
        grid.expand_volume(3);
        hydration_atoms = grid.hydrate();

        // double width = 10; // what width to use? 10 is too large, but with smaller values our grid becomes incredibly large
        // auto[corner, bins] = generate_grid(width); // corner is the lower corner of our grid, and bins the number of bins in each dimension
        // cout << format("bins: (%1%, %2%, %3%)") % bins[0] % bins[1] % bins[2] << endl;

        // Jans approach: generate grid from -250 to 250 in all dimensions, with grid size 1Å. 
        // Move protein to center so its more likely to be covered by the grid. 
        // Determine which grid points are filled by our atoms. We can do this in a smart way by converting each coordinate to an index instead of searching through the entire grid. 
        // Expand each grid point into a proper size for its atom. Jan uses sqrt(6) in each direction. 
        // Calculate the volume by iterating through the grid and counting how many entries are filled. Can probably be done smarter by keeping track of it along the way.
        //      The volume may be scaled. Jan multiplied by 27 for some reason. 
        // Determine the surface atoms of the protein. Jan did this by checking 3Å in each direction. If no other atom occupies the spot, it is a candidate for a hydration atom. 
        // Place the hydration atoms. Only keep every third. 
    }

    vector<shared_ptr<Atom>>* get_protein_atoms() {
        return &protein_atoms;
    }

    vector<shared_ptr<Atom>>* get_hydration_atoms() {
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
        auto update = [&low, &high] (vector<shared_ptr<Atom>>* atoms) {
            for (auto const& a : *atoms) {
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
        auto weighted_sum = [&cm, &M] (vector<shared_ptr<Atom>>* atoms) {
            for (auto const& a : *atoms) {
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

private:
    vector<shared_ptr<Atom>> protein_atoms; // atoms of the protein itself
    vector<shared_ptr<Atom>> hydration_atoms; // hydration layer
    shared_ptr<File> file; 

    /** Move the entire protein by a vector
     * @param v the translation vector
     */
    void translate(const TVector3 v) {
        auto move = [&v] (vector<shared_ptr<Atom>>* atoms) {
            for (auto const& a : *atoms) {
                a->translate(v);
            }
        };
        move(&protein_atoms);
        move(&hydration_atoms);
    }
};