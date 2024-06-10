/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/files/Folder.h>

#include <filesystem>

using namespace md;

Folder::Folder() = default;

Folder::Folder(const std::string& path) : path(path) {
    if (path != "") {std::filesystem::create_directories(path);}
}

Folder::operator std::string() const {
    return path;
}

std::string Folder::operator+(const std::string& s) const {
    return path + s;
}