/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/Folder.h>
#include <utility/Exceptions.h>

#include <filesystem>

using namespace io;


Folder::Folder(const char* path) : Folder(std::string(path)) {}
Folder::Folder(const std::string& path) {
    *this = path;
}

void Folder::operator=(const std::string& path) {
    if (path.empty()) {dir = "."; return;}
    if (path.back() == '/') {
        dir = path.substr(0, path.size() - 1);
    } else {
        dir = path;
    }
}

Folder& Folder::operator+(const std::string& str) noexcept {
    auto s = str.back() == '/' ? str.substr(0, str.size() - 1) : str;
    dir += s.front() == '/' ? std::move(s) : "/" + std::move(s);
    return *this;
}

Folder& Folder::operator+(const Folder& folder) noexcept {
    return *this + std::string(folder);
}

Folder::operator std::string() const {return dir;}

std::string Folder::path() const {return dir;}

bool Folder::exists() const noexcept {
    if (dir.empty()) {return false;}
    return std::filesystem::exists(dir);
}

std::vector<std::string> Folder::files() const {
    throw except::not_implemented("io::Folder::files");
}

std::vector<std::string> Folder::directories() const {
    throw except::not_implemented("io::Folder::directories");
}

void Folder::create() const {
    if (exists()) {return;}
    std::filesystem::create_directories(dir);
}

std::string operator+(const char* str, const io::Folder& folder) {
    return std::string(str) + folder;
}

std::string operator+(const std::string& str, const io::Folder& folder) {
    auto s = str.back() == '/' ? str.substr(0, str.size() - 1) : str;
    if (folder.path().front() == '/') {return std::move(s) + folder.path();}
    else {return std::string(std::move(s)) + "/" + folder.path();}
}