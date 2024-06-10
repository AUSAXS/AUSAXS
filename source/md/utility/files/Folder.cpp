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