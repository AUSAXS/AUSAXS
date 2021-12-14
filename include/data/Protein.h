// includes
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

// my own stuff
#include "data/Body.h"
#include "Tools.h"
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "io/File.h"

using std::vector, std::string;

class Protein {
public: 
    Protein(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms);
    Protein(const string& input);

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

private:
    vector<Body> bodies;
    std::unique_ptr<File> file;

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
    void calc_distances();

    /**
     * @brief Subtract the charge of the displaced water molecules from the effective charge of the protein atoms. 
     */
    void update_effective_charge();
};