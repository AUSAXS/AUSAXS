// includes
#include <string>
#include <vector>
#include <utility>

// my own stuff
#include "data/Body.h"
#include "Tools.h"
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "data/StateManager.h"
#include "data/PartialHistogramManager.h"

using std::vector, std::string;

class Protein {
public: 
    Protein() {}

    /**
     * @brief Create a new protein based on vectors of atoms.
     */
    Protein(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms);

    /**
     * @brief Create a new protein from a list of input sources.
     * @param input a list of paths to the input files. File extensions can be mixed. 
     */
    Protein(const vector<string>& input);

    /**
     * @brief Create a new protein from a single input source. 
     * @param input the path to the input file. 
     */
    Protein(const string& input);

    /**
     * @brief Get the distances between each atom.
     */
    shared_ptr<ScatteringHistogram> get_histogram();

    /** 
     * @brief Writes this body to disk.
     * @param path path to the destination. 
     */
    void save(string path);

    /** 
     * @brief Use an algorithm to generate a new hydration layer for this body. Note that the previous one will be deleted.
     */
    void generate_new_hydration();

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
     * @brief Calculate the center-mass coordinates for the protein.
     * @return The center-mass (x, y, z) coordinates. 
     */
    Vector3 get_cm() const;

    /**
     * @brief Calculate the total mass of this protein in Daltons.
     */
    double get_mass() const;

    /**
     * @brief Get a copy of all protein atoms from the underlying bodies.
     */
    vector<Atom> get_protein_atoms() const;

    /**
     * @brief Get a copy of the hydration atoms. Use the member variable for reference access. 
     */
    vector<Hetatom> get_hydration_atoms() const;

    vector<Hetatom> hydration_atoms; // stores the hydration atoms from the generated hydration layer
private:
    vector<Body> bodies; // the constituent bodies
    shared_ptr<Grid> grid = nullptr; // the grid representation of this body
    unique_ptr<PartialHistogramManager> phm = nullptr;
    shared_ptr<ScatteringHistogram> histogram = nullptr; // an object representing the distances between atoms

    /** 
     * @brief Move the entire protein by a vector.
     * @param v the translation vector.
     */
    void translate(const Vector3& v);

    /**
     * @brief Create a grid and fill it with the atoms of this protein. 
     */
    void create_grid();

    /** 
     * @brief Calculate the distances between each pair of atoms. 
     */
    void calc_histogram();

    /**
     * @brief Subtract the charge of the displaced water molecules from the effective charge of the protein atoms. 
     */
    void update_effective_charge();
};