// includes
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

// my own stuff
#include "data/Record.h"
#include "Tools.h"
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "io/File.h"

using std::vector, std::string;

class Body {
public: 
    Body(const vector<Atom>& protein_atoms, const vector<Hetatom>& hydration_atoms);
    Body(const string& input);

    void translate(const Vector3& v);

    /**
     * @brief Rotatate the Body @a rad radians about the axis @a axis. 
     * @param axis the rotation axis. 
     * @param rad the amount to rotate in radians. 
     */
    void rotate(const Vector3& axis, const double& rad);

    /**
     * @brief Euler angle rotation. 
     * @param alpha radians to rotate about the z-axis.
     * @param beta radians to rotate about the y-axis. 
     * @param gamma radians to rotate about the x-axis. 
     */
    void rotate(const double& alpha, const double& beta, const double& gamma);

private:
    vector<Atom>& protein_atoms;
    vector<Hetatom>& hydration_atoms;
    std::unique_ptr<File> file;
};