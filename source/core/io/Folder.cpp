// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/Folder.h>

#include <io/File.h>

#include <filesystem>

using namespace ausaxs;
using namespace ausaxs::io;

Folder::Folder(std::string_view path) {*this = path;}

Folder& Folder::operator=(std::string_view path) {
    if (path.empty()) {dir = "."; return *this;}
    auto is_separator = [] (char c) {return c == std::filesystem::path::preferred_separator;};
    while (1 < path.size() && is_separator(path.back())) {path.remove_suffix(1);} // keep the root folder intact
    dir = path;
    return *this;
}

Folder::operator std::string() const {return dir;}

std::string Folder::path() const {return dir;}

std::string Folder::str() const {return path();}

bool Folder::empty() const noexcept {return dir.empty();}

bool Folder::exists() const noexcept {
    if (dir.empty()) {return false;}
    std::error_code ec;
    return std::filesystem::is_directory(dir, ec);
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

void Folder::remove() const {
    if (dir.empty()) {return;}
    std::filesystem::remove_all(dir);
}

io::Folder Folder::move(const io::Folder& destination) const {
    std::filesystem::rename(dir, destination.path());
    return destination;
}

template Folder::Folder(const char* const&);
template Folder::Folder(const std::string&);
template Folder::Folder(const std::string_view&);
template Folder& Folder::operator=(const std::string&);
template Folder& Folder::operator=(const std::string_view&);