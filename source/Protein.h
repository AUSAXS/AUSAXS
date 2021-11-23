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
#include "data/Distances.h"

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
    void save(string path) const;

    /**
     * @brief Get the distances between each atom.
     */
    shared_ptr<Distances> get_distances();

    /**
     * @brief Calculate the intensity based on the Debye scattering equation
     */
    vector<double> debye_scattering_intensity() const;

    /**
     * @brief Calculate the intensity based on the Debye scattering equation
     * @param axes the axes for p in the format {bins, xmin, xmax}
     * @param p the binned distance histogram to use
     */
    vector<double> debye_scattering_intensity(vector<int> axes, vector<double>& p) const;

    /** 
     * @brief Use an algorithm to generate a new hydration layer for this protein. Note that the previous one will be deleted.
     */
    void generate_new_hydration();

    /**
     * @brief Get a pointer to the protein atoms.
     */
    vector<shared_ptr<Atom>>* get_protein_atoms() {return &protein_atoms;}

    /**
     * @brief Get a pointer to the hydration atoms.
     */
    vector<shared_ptr<Hetatom>>* get_hydration_atoms() {return &hydration_atoms;}

    /** Calculate the center-mass coordinates for the protein.
     * @return The center-mass (x, y, z) coordinates. 
     */
    TVector3 get_cm() const;

    /**
     * @brief Calculate the volume of this protein based on its constituent amino acids
     */
    double get_volume() const;

    /**
     * @brief Get the grid representation of this Protein. 
     */
    shared_ptr<Grid> get_grid() const {
        if (grid == nullptr) {
            print_err("Error in Protein::get_grid: Grid has not been instantiated!"); 
            exit(1);
        }
        return grid;
    }

    /**
     * @brief Generate a PDB file showing the filled grid volume.
     */
    void generate_volume_file(string path);

private:
    vector<shared_ptr<Atom>> protein_atoms; // atoms of the protein itself
    vector<shared_ptr<Hetatom>> hydration_atoms; // hydration layer
    shared_ptr<File> file; // the file backing this protein
    shared_ptr<Grid> grid = nullptr; // the grid representation of this protein
    shared_ptr<Distances> distances; // an object representing the distances between atoms

    /** Move the entire protein by a vector.
     * @param v the translation vector
     */
    void translate(const TVector3 v);

    /** 
     * @brief Calculate the distances between each pair of atoms. 
     */
    void calc_distances();
};