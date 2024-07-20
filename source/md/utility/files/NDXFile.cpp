/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/files/NDXFile.h>

#include <fstream>

using namespace md;

void NDXFile::append_file(const NDXFile& other) {
    std::ifstream ifs(other);
    std::ofstream ofs(*this, std::ios::app);
    ofs << ifs.rdbuf();
}

bool NDXFile::contains(const std::string& group) const {
    std::ifstream ifs(*this);
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) {continue;}
        if (line[0] != '[') {continue;}
        auto end = line.find(']');
        if (end == std::string::npos) {throw std::runtime_error("Invalid index file \"" + path() + "\": line \"" + line + "\" does not contain a closing bracket.");}
        if (line[1] != ' ' || line[end-1] != ' ') {throw std::runtime_error("Invalid index file\"" + path() + "\": line \"" + line + "\" does not contain a space before and after the group name.");}
        line = line.substr(2, end-3);
        if (line == group) {
            return true;
        }
    }
    return false;
}