#include <io/File.h>
#include <utility/Exceptions.h>

#include <filesystem>
#include <fstream>

using namespace io;

File::File(const io::Folder& folder, std::string_view name, std::string_view extension) : dir(folder), name(name), ext(extension) {}

File::File(const char* path) : File(std::string(path)) {}

File::File(const std::string& path) {
    if (path.empty()) {return;}
    *this = path;
}

std::tuple<std::string, std::string, std::string> File::split(const std::string& path) {
    auto p = std::filesystem::path(path);
    return std::make_tuple(p.parent_path(), p.stem(), p.extension());
}

void File::operator=(const std::string& path) {
    auto [dir, file, ext] = split(path);
    this->dir = std::move(dir);
    this->name = std::move(file);
    this->ext = std::move(ext);
}

File::operator std::string() const {
    return dir.path() + "/" + name + ext;
}

File& File::replace_extension(const std::string& extension) noexcept {
    if (extension.front() == '.') {
        ext = extension.substr(1);
    } else {
        ext = extension;
    }
    return *this;
}

File& File::append(const std::string& name) noexcept {
    this->name += name;
    return *this;
}

io::File File::append(const std::string& name) const noexcept {
    auto file = *this;
    file.append(name);
    return file;
}

Folder& File::directory() noexcept {return dir;}

const Folder& File::directory() const noexcept {return dir;}

std::string& File::extension() noexcept {return ext;}

const std::string& File::extension() const noexcept {return ext;}

std::string File::stem() const noexcept {return name;}

void File::create(const std::string& contents) const {
    if (!dir.exists()) {dir.create();}
    std::ofstream file(*this);

    if (!file.is_open()) {throw except::io_error("File::create: Could not create file " + path());}
    file << contents;
}

void File::remove() const {
    if (exists()) {std::filesystem::remove(path());}
}

std::string File::path() const {return std::string(*this);}

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