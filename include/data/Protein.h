#pragma once

// forwards declaration
class ScatteringHistogram;

// includes
#include <vector>
#include <map>
#include <utility>

// my own includes
#include "data/Atom.h"
#include "hydrate/Grid.h"
#include "io/File.h"
#include "constants.h"
#include "ScatteringHistogram.h"

using std::vector, std::string, std::unique_ptr;
using namespace ROOT;

class Protein {
public:
    /** Creates a new protein from the input .pdb or .xml file. 
     * @param path path to the input file. 
     */
    Protein(const string& path) : file(std::make_shared<File>(path)), protein_atoms(file->protein_atoms), hydration_atoms(file->hydration_atoms) {}

    /** Writes this protein to disk.
     * @param path path to the destination. 
     */
    void save(string path);

    /**
     * @brief Get the distances between each atom.
     */
    shared_ptr<ScatteringHistogram> get_distances();

    /** 
     * @brief Use an algorithm to generate a new hydration layer for this protein. Note that the previous one will be deleted.
     */
    void generate_new_hydration();

    /**
     * @brief Get a pointer to the protein atoms.
     */
    const vector<Atom>& get_protein_atoms() {return protein_atoms;}

    /**
     * @brief Get a pointer to the hydration atoms.
     */
    const vector<Hetatom>& get_hydration_atoms() {return hydration_atoms;}

    /** Calculate the center-mass coordinates for the protein.
     * @return The center-mass (x, y, z) coordinates. 
     */
    Vector3 get_cm() const;

    /**
     * @brief Calculate the volume of this protein based on its constituent amino acids
     */
    double get_volume_acids() const;

    /**
     * @brief Calculate the volume of this protein based on the number of grid bins it spans
     */
    double get_volume_grid();

    /**
     * @brief Calculate the volume of this protein based on the number of C-alpha atoms
     */
    double get_volume_calpha() const;

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
     * @brief Generate a PDB file at @p path showing the filled grid volume.
     */
    void generate_volume_file(string path);

    /**
     * @brief Calculate the total mass of this protein.
     */
    double get_mass() const;

private:
    shared_ptr<File> file = nullptr; // the file backing this protein
    vector<Atom>& protein_atoms; // atoms of the protein itself
    vector<Hetatom>& hydration_atoms; // hydration layer
    shared_ptr<Grid> grid = nullptr; // the grid representation of this protein
    shared_ptr<ScatteringHistogram> distances = nullptr; // an object representing the distances between atoms

    /**
     * @brief Create a grid and fill it with the atoms of this protein. 
     */
    void create_grid();

    /** Move the entire protein by a vector.
     * @param v the translation vector
     */
    void translate(const Vector3& v);

    /** 
     * @brief Calculate the distances between each pair of atoms. 
     */
    void calc_distances();

    /**
     * @brief Subtract the charge of the displaced water molecules from the effective charge of the protein atoms. 
     */
    void update_effective_charge();
};