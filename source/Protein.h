#pragma once

// includes
#include <vector>
#include <map>
#include <utility>

// ROOT
#include <TVector3.h>

// my own includes
#include "data/Atom.h"
#include "hydrate/Grid.h"
#include "data/File.h"
#include "data/properties.h"

using std::vector, std::string, std::unique_ptr;
using namespace ROOT;

class Protein {
public:
    /** Creates a new protein from the input .pdb or .xml file. 
     * @param path path to the input file. 
     */
    Protein(string path);

    /** Writes this protein to disk.
     * @param path path to the destination. 
     */
    void save(string path);

    /** Calculate the distances between each pair of atoms. 
     * @return A pair where the first entry is a vector of all internal distances between the protein atoms, while the second entry is all internal
     * distances between hydration atoms plus distances between hydration and protein atoms. 
     */
    std::pair<vector<double>, vector<double>> calc_distances();

    /** 
     * @brief Use an algorithm to generate a new hydration layer for this protein. Note that the previous one will be deleted.
     * @param reduce the factor to reduce the output number of water molecules by. 
     * @param width the distance between each grid point
     */
    void generate_new_hydration(int reduce, double width);

    /**
     * @brief Get a pointer to the protein atoms.
     */
    vector<shared_ptr<Atom>>* get_protein_atoms() {return &protein_atoms;}

    /**
     * @brief Get a pointer to the hydration atoms.
     */
    vector<shared_ptr<Hetatom>>* get_hydration_atoms() {return &hydration_atoms;}

    /** Generate a 3D grid containing all atoms.
     * @param width the bin width
     * @return A pair (corner point, bins) where the first is a bottom corner of the grid, while the second is the number of bins in each dimension.
     */
    std::pair<TVector3, vector<int>> generate_grid(const double width);

    /** Calculate the center-mass coordinates for the protein.
     * @return The center-mass (x, y, z) coordinates. 
     */
    TVector3 get_cm();

    /**
     * @brief Calculate the volume of this protein based on its constituent amino acids
     */
    double get_volume();

private:
    vector<shared_ptr<Atom>> protein_atoms; // atoms of the protein itself
    vector<shared_ptr<Hetatom>> hydration_atoms; // hydration layer
    shared_ptr<File> file; 

    /** Move the entire protein by a vector
     * @param v the translation vector
     */
    void translate(const TVector3 v);
};