#pragma once

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
    /**
     * @brief Default constructor. 
     */
    Protein() {}

    /**
     * @brief Constructor.
     * 
     * Create a new protein based on a set of bodies.
     * 
     * @param bodies The constituent bodies of this protein. 
     * @param hydration_atoms The hydration layer. 
     */
    Protein(const vector<Body>& bodies, const vector<Hetatom>& hydration_atoms = {});

    /**
     * @brief Constructor.
     * 
     * Create a new protein based on a set of atoms. 
     * This will only create a single constituent body. 
     * 
     * @param protein_atoms The constituent atoms of this protein. 
     * @param hydration_atoms The hydration layer. 
     */
    Protein(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms = {});

    /**
     * @brief Constructor. 
     * 
     * Create a new protein based on a set of atom vectors. Each vector defines a constituent body. 
     * 
     * @param protein_atoms The constituent atoms of each body. 
     * @param hydration_atoms The hydration layer. 
     */
    Protein(const vector<vector<Atom>>& protein_atoms, const vector<Hetatom>& hydration_atoms = {});

    /**
     * @brief Constructor. 
     * 
     * Create a new protein based on a list of input file paths. 
     * 
     * @param input A list of paths to the input files. File extensions can be mixed. 
     */
    Protein(const vector<string>& input);

    /**
     * @brief Constructor.
     * 
     * Create a new protein based on a single input file path. 
     * 
     * @param input Path to the input file. 
     */
    Protein(const string& input);

    /**
     * @brief Get the distances between each atom.
     */
    shared_ptr<ScatteringHistogram> get_histogram();

    /**
     * @brief Get the total distance histogram only. 
     *        This is a slightly faster alternative to get_histogram() when only the total histogram is needed. 
     */
    Histogram get_total_histogram() const;

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
     * @brief Get the grid representation of this body. 
     */
    shared_ptr<Grid> get_grid() const;

    /**
     * @brief Create a binding point between two bodies.
     *        This binding will be a constraint for rigid-body optimization. 
     */
    void bind();

    /**
     * @brief Center this protein on origo. 
     */
    void center();

    /**
     * @brief Get a copy of all protein atoms from the underlying bodies.
     */
    vector<Atom> get_protein_atoms() const;

    /**
     * @brief Get a copy of the hydration atoms. Use the member variable for reference access. 
     */
    vector<Hetatom> get_hydration_atoms() const;

    vector<Hetatom> hydration_atoms; // Stores the hydration atoms from the generated hydration layer
    vector<Body> bodies; // The constituent bodies
    bool updated_charge = false; // True if the effective charge of each atom has been updated to reflect the volume they occupy, false otherwise
    bool centered = false; // True if this object is centered, false otherwise. 
  private:
    shared_ptr<Grid> grid = nullptr; // The grid representation of this body
    unique_ptr<PartialHistogramManager> phm = nullptr;
    shared_ptr<ScatteringHistogram> histogram = nullptr; // An object representing the distances between atoms

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