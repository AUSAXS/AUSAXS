#include <md/utility/files/File.h>
#include <md/utility/Exceptions.h>

#include <filesystem>
#include <fstream>

using namespace md::detail;

File::File() = default;
File::File(const std::string& path) : path(path) {}
File::File(const std::string& path, const std::string& ext) : path(path), ext(ext) {
    validate(path);
}

void File::validate(const std::string& path) {
    if (path.substr(path.find_last_of(".") + 1) != ext) {
        throw except::invalid_format("gmx::File: \"" + path + "\" does not have the correct extension (" + ext + ").");
    }
}

File::operator std::string() const {
    return path;
}

bool File::empty() const {
    return path.empty();
}

std::string File::parent_path() const {
    if (empty()) {throw except::io_error("File::parent_path: Cannot get parent path of empty file.");}
    return std::filesystem::path(path).parent_path().string();
}

void File::remove() const {
    if (!exists()) {throw except::io_error("File::remove: Cannot remove file \"" + path + "\" because it does not exist.");}
    std::filesystem::remove(path);
}

void File::create(const std::string& content) const {
    if (exists()) {throw except::io_error("File::create: Cannot create file \"" + path + "\" because it does not exist.");}
    std::ofstream file(path);
    if (!file.is_open()) {
        throw except::io_error("gmx::File: Could not create file \"" + path + "\".");
    }
    file << content;
}

std::string File::move(const Folder& folder) {
    if (!std::filesystem::exists(path)) {throw except::io_error("File::move: Cannot move file \"" + path + "\" because it does not exist.");}

    std::filesystem::rename(path, folder + filename());
    path = folder + filename();
    return folder + filename();
}

std::string File::copy(const Folder& folder) const {
    if (!std::filesystem::exists(path)) {throw except::io_error("File::copy: Cannot copy file \"" + path + "\" because it does not exist.");}

    std::string newfile = folder;
    if (newfile.back() != '/') {newfile += "/";}
    newfile += filename();
    std::filesystem::copy_file(path, newfile, std::filesystem::copy_options::overwrite_existing);
    return newfile;
}

File& File::rename(const std::string& newname) {
    if (!std::filesystem::exists(path)) {throw except::io_error("File::rename: Cannot rename file \"" + path + "\" because it does not exist.");}

    std::filesystem::rename(path, parent_path() + "/" + newname);
    path = parent_path() + "/" + newname;
    return *this;
}

std::string File::filename() const {
    if (empty()) {throw except::io_error("File::filename: Cannot get filename of empty path.");}
    return std::filesystem::path(path).filename().string();
}

bool File::exists() const {
    if (empty()) {throw except::io_error("File::exists: Cannot check if empty path exists.");}
    return std::filesystem::exists(path);
}

std::string File::relative_path(const File& other) const {
    if (empty()) {throw except::io_error("File::relative_path: Cannot get relative path of empty path.");}
    return relative_path(other.path);
}

std::string File::relative_path(const std::string& s) const {
    if (empty()) {throw except::io_error("File::relative_path: Cannot get relative path of empty path.");}
    return std::filesystem::relative(s, parent_path()).string();
}

std::string File::absolute() const {
    if (empty()) {throw except::io_error("File::absolute: Cannot get absolute path of empty path.");}
    return std::filesystem::absolute(path).string();
}

std::string File::stem() const {
    if (empty()) {throw except::io_error("File::stem: Cannot get stem of empty path.");}
    return std::filesystem::path(path).stem().string();
}