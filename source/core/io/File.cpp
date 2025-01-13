/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/File.h>
#include <utility/Exceptions.h>

#include <filesystem>
#include <fstream>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::io;

File::File(const io::Folder& folder, std::string_view name, std::string_view extension) : dir(folder), name(name), ext(extension) {}
File::File(std::string_view name, std::string_view extension) : File(Folder(), name, extension) {}

File::File(std::string_view path) {
    auto[dir, file, ext] = split(path);
    this->dir = std::move(dir);
    this->name = std::move(file);
    this->ext = std::move(ext);
}

std::tuple<std::string, std::string, std::string> File::split(std::string_view path) {
    auto p = std::filesystem::path(path);
    return std::make_tuple(p.parent_path().string(), p.stem().string(), p.extension().string());
}

File::operator std::string() const {
    return dir.path() + "/" + name + ext;
}

File& File::replace_extension(std::string_view extension) noexcept {
    if (extension.front() == '.') {
        ext = extension.substr(1);
    } else {
        ext = extension;
    }
    return *this;
}

File& File::append(std::string_view name) noexcept {
    this->name += name;
    return *this;
}

io::File File::append(std::string_view name) const noexcept {
    auto file = *this;
    file.append(name);
    return file;
}

Folder& File::directory() noexcept {return dir;}

const Folder& File::directory() const noexcept {return dir;}

std::string& File::extension() noexcept {return ext;}

const std::string& File::extension() const noexcept {return ext;}

std::string File::stem() const noexcept {return name;}

std::string& File::stem() noexcept {return name;}

void File::create(std::string_view contents) const {
    if (!dir.exists()) {dir.create();}
    std::ofstream file(*this);

    if (!file.is_open()) {throw except::io_error("File::create: Could not create file " + path());}
    file << contents;
}

io::File File::rename(std::string_view name) const {
    io::File new_file(dir, name, ext);
    std::filesystem::rename(path(), new_file.path());
    return new_file;
}

io::File File::move(const io::Folder& folder) const {
    if (folder == dir) {return *this;}
    if (!folder.exists()) {folder.create();}
    io::File new_file(folder, name, ext);
    if (new_file.exists()) {new_file.remove();}
    std::filesystem::rename(path(), new_file.path());
    return new_file;
}

io::File File::copy(const io::Folder& folder) const {
    if (folder == dir) {return *this;}
    if (!folder.exists()) {folder.create();}
    io::File new_file(folder, name, ext);
    if (new_file.exists()) {new_file.remove();}
    std::filesystem::copy(path(), new_file.path());
    return new_file;
}

void File::remove() const {
    if (exists()) {std::filesystem::remove(path());}
}

std::string File::path() const {return std::string(*this);}
std::string File::str() const {return path();}

std::string File::absolute_path() const {return std::filesystem::absolute(path()).string();}

std::string File::relative_path(const File& other) const {
    if (name.empty() || other.name.empty()) {throw except::io_error("File::relative_path: Cannot get relative path of empty path.");}
    return std::filesystem::relative(other.absolute_path(), directory().path()).string();
}

std::string File::filename() const noexcept {
    return name + ext;
}

bool File::empty() const noexcept {
    return name.empty() && ext.empty();
}

bool File::exists() const noexcept {
    return std::filesystem::is_regular_file(path());
}

std::ostream& operator<<(std::ostream& os, const io::File& file) {
    os << file.path();
    return os;
}