/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/Protein.h>
#include <md/utility/Exceptions.h>

#include <fstream>
#include <cmath>

using namespace md;

std::vector<Vector3> Protein::parse(const std::string& filename) {
    std::ifstream input(filename);
    if (!input.is_open()) {throw except::io_error("PDBParser::parse: Could not open file \"" + filename + "\".");}

    std::vector<Vector3> result;
    std::string line;
    while(std::getline(input, line)) {
        if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
            result.push_back({std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)), std::stod(line.substr(46, 8))});
        }
    }
    return result;
};

double Protein::maximum_distance() const {
    double result = 0;
    for (auto i = atoms.begin(); i != atoms.end(); ++i) {
        for (auto j = i + 1; j != atoms.end(); ++j) {
            double distance = std::sqrt(std::pow(i->x - j->x, 2) + std::pow(i->y - j->y, 2) + std::pow(i->z - j->z, 2));
            if (distance > result) {result = distance;}
        }
    }
    return result;
}