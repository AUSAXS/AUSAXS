#include <io/File.h>

#include <filesystem>
#include <fstream>

using namespace io;

File::File(const std::string& path) {
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
    return dir.path() + "/" + name + "." + ext;
}

void File::replace_extension(const std::string& extension) noexcept {
    if (extension.front() == '.') {
        ext = extension.substr(1);
    } else {
        ext = extension;
    }
}

void File::append(const std::string& name) noexcept {
    this->name += name;
}

Folder& File::directory() noexcept {return dir;}

const Folder& File::directory() const noexcept {return dir;}

std::string& File::extension() noexcept {return ext;}

const std::string& File::extension() const noexcept {return ext;}

std::string File::stem(const std::string& path) {return name;}

void File::create() const {
    if (!dir.exists()) {dir.create();}
    std::ofstream file(*this);
}

std::string File::path() const {return std::string(*this);}

std::string operator+(const char* str, const io::File& file) {
    return std::string(str) + std::string(file);
}

bool File::exists() const noexcept {
    std::filesystem::exists(path());
}