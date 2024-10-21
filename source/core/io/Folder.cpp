/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/Folder.h>
#include <io/File.h>
#include <utility/Exceptions.h>

#include <filesystem>

using namespace io;

Folder::Folder(std::string_view path) {*this = path;}

void Folder::operator=(std::string_view path) {
    if (path.empty()) {dir = "."; return;}
    if (path.back() == '/') {
        dir = path.substr(0, path.size() - 1);
    } else {
        dir = path;
    }
}

Folder::operator std::string() const {return dir;}

std::string Folder::path() const {return dir;}

bool Folder::empty() const noexcept {return dir.empty();}

bool Folder::exists() const noexcept {
    if (dir.empty()) {return false;}
    return std::filesystem::is_directory(dir);
}

std::vector<io::File> Folder::files() const {
    std::vector<io::File> files;
    for (const auto& entry : std::filesystem::directory_iterator(path())) {
        if (entry.is_regular_file()) {files.emplace_back(entry.path().string());}
    }
    return files;
}

std::vector<io::Folder> Folder::directories() const {
    std::vector<io::Folder> dirs;
    for (const auto& entry : std::filesystem::directory_iterator(path())) {
        if (entry.is_directory()) {dirs.emplace_back(entry.path().string());}
    }
    return dirs;
}

void Folder::create() const {
    if (exists()) {return;}
    std::filesystem::create_directories(dir);
}

std::string operator+(const io::Folder& folder, std::string_view str) {
    return folder.path() + "/" + std::string(str);
}

std::string operator+(std::string_view str, const io::Folder& folder) {
    auto s = str.back() == '/' ? str.substr(0, str.size() - 1) : str;
    return std::string(s) + "/" + folder.path();
}

template Folder::Folder(const char* const&);
template Folder::Folder(const std::string&);
template Folder::Folder(const std::string_view&);
template void Folder::operator=(const std::string&);
template void Folder::operator=(const std::string_view&);