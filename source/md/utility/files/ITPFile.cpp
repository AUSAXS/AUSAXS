/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/files/ITPFile.h>
#include <md/utility/Exceptions.h>
#include <utility/StringUtils.h>

#include <fstream>

using namespace md;

unsigned int ITPFile::size() const {
    std::ifstream in(path());
    std::string line, last_atom_line;
    bool in_atom_section = false;
    unsigned int count = 0;
    while (std::getline(in, line)) {
        if (line[0] == ';' || line.empty()) {
            continue;
        }

        if (line.find("[ atoms ]") == 0) {
            in_atom_section = true;
            continue;
        } else if (in_atom_section && (line[0] == '[' || line[0] == '#')) {
            in_atom_section = false;
        }

        if (in_atom_section) {
            count++;
            last_atom_line = line;
        }
    }

    // check that the count & last atom line agrees
    if (count != 0) {
        auto tokens = utility::split(last_atom_line, " \t");
        try {
            if (std::stoi(tokens[0]) != static_cast<int>(count)) {
                throw except::invalid_argument("ITPFile::size: inconsistent atom count in file \"" + path() + "\". Expected " + std::to_string(count) + " but found " + tokens[0] + ".");
            }
        } 
        catch (const std::invalid_argument& e) {
            throw except::invalid_argument("ITPFile::size: could not parse atom index from line \"" + last_atom_line + "\" in file \"" + path() + "\".");
        }
    }

    if (1e6 < count) {
        throw except::invalid_argument("ITPFile::size: atom count in file \"" + path() + "\" seems unreasonably large (" + std::to_string(count) + "). Contact the developer if this is an error.");
    }

    return count;
}

std::vector<ITPFile> ITPFile::split_restraints(const std::vector<ITPFile>& topologies) const {
    std::ifstream in(path());
    if (!in.is_open()) {throw except::io_error("ITPFile::split_restraints: could not open file \"" + path() + "\".");}
    if (topologies.size() == 0) {throw except::invalid_argument("ITPFile::split_restraints: topologies must not be empty.");}
    for (const auto& topol : topologies) {
        if (topol.path().find("topol") == std::string::npos) {
            throw except::invalid_argument("ITPFile::split_restraints: topologies must be the original include topologies, named topol_*.itp.");
        }
    }

    std::vector<ITPFile> files;
    std::vector<unsigned int> split_indices = {0};

    unsigned int sum = 0;
    for (const auto& topol : topologies) {
        sum += topol.size();
        split_indices.push_back(sum);
    }

    std::vector<std::vector<std::string>> file_contents(topologies.size());
    for (unsigned int i = 0; i < topologies.size(); i++) {
        file_contents[i].push_back(
            "[ position_restraints ]\n"
            ";  i funct       fcx        fcy        fcz"
        );
    }

    unsigned int current_topology = 1;

    std::string line;
    while (std::getline(in, line)) {
        if (line[0] == ';' || line.empty()) {
            continue;
        }

        if (line[0] == '[') {
            continue;
        }

        auto tokens = utility::split(line, " \t");
        unsigned int index;
        try {index = std::stoi(tokens[0]);} 
        catch (const std::invalid_argument& e) {throw except::invalid_argument("ITPFile::split_restraints: could not parse atom index from line \"" + line + "\" in file \"" + path() + "\".");}
        if (index >= split_indices[current_topology]) {
            current_topology++;
            if (current_topology == topologies.size()) {
                throw except::invalid_argument("ITPFile::split_restraints: atom index " + std::to_string(index) + " exceeds the number of atoms in the topologies (" + std::to_string(split_indices.back()) + ").");
            }
        }

        index -= split_indices[current_topology-1];
        tokens[0] = std::to_string(index);
        line = utility::join(tokens, "\t");

        file_contents[current_topology-1].push_back(line);
    }

    for (unsigned int i = 0; i < topologies.size(); i++) {
        // replace "topol" from filename with "backbone" 
        std::string filename = topologies[i];
        filename.replace(filename.find("topol"), 5, "backbone");
        std::ofstream out(filename);
        for (const auto& line : file_contents[i]) {
            out << line << "\n";
        }
        files.push_back(ITPFile(filename));
    }

    return files;
}