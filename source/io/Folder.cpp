#include <io/Folder.h>

#include <filesystem>

using namespace io;

Folder::Folder(const std::string& path) : dir(path) {
    *this = path;
}

void Folder::operator=(const std::string& path) {
    bool slash_b = path.front() == '/';
    bool slash_f = path.back() == '/';
    if (slash_b && slash_f) {
        dir = path.substr(1, path.size() - 2);
    } else if (slash_b) {
        dir = path.substr(1);
    } else if (slash_f) {
        dir = path.substr(0, path.size() - 1);
    } else {
        dir = path;
    }
}

Folder& Folder::operator+(const Folder& folder) noexcept {
    dir += "/" + folder.dir;
    return *this;
}

Folder::operator std::string() const {return dir;}

std::string Folder::path() const {return dir;}

bool Folder::exists() const noexcept {
    if (dir.empty()) {return false;}
    return std::filesystem::exists(dir);
}

std::vector<std::string> Folder::files() const {
    std::vector<std::string> files;
    return files;
}

std::vector<std::string> Folder::directories() const {
    std::vector<std::string> directories;
    return directories;
}

void Folder::create() const {
    if (exists()) {return;}
    std::filesystem::create_directories(dir);
}