/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/File.h>
#include <utility/Exceptions.h>

#include <filesystem>
#include <fstream>

using namespace io;

File::File(const io::Folder& folder, std::string_view name, std::string_view extension) : dir(folder), name(name), ext(extension) {}

File::File(std::string_view path) {
    if (path.empty()) {return;}

    auto [dir, file, ext] = split(path);
    this->dir = std::move(dir);
    this->name = std::move(file);
    this->ext = std::move(ext);
}

File::File(const char* path) : File(std::string_view(path)) {}

File::File(const std::string& path) : File(std::string_view(path)) {}

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

void File::create(std::string_view contents) const {
    if (!dir.exists()) {dir.create();}
    std::ofstream file(*this);

    if (!file.is_open()) {throw except::io_error("File::create: Could not create file " + path());}
    file << contents;
}

void File::remove() const {
    if (exists()) {std::filesystem::remove(path());}
}

std::string File::path() const {return std::string(*this);}

std::string File::absolute_path() const {return std::filesystem::absolute(path()).string();}

bool File::exists() const noexcept {
    return std::filesystem::exists(path());
}

std::string operator+(const char* str, const io::File& file) {
    return std::string(str) + std::string(file);
}

std::string operator+(const io::File& file, const char* str) {
    return std::string(file) + std::string(str);
}

std::string operator+(const std::string& str, const io::File& file) {
    return str + std::string(file);
}

std::string operator+(const io::File& file, const std::string& str) {
    return std::string(file) + str;
}

std::istringstream &operator>>(std::istringstream& in, io::File& val) {
    std::string v;
    in >> v;
    val = v;
    return in;
}

std::ostream& operator<<(std::ostream& os, const io::File& file) {
    os << file.path();
    return os;
}