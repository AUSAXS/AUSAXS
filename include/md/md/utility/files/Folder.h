#pragma once

#include <string>

namespace gmx {
    struct Folder {
        Folder();
        Folder(const std::string& path);
        operator std::string() const;
        std::string operator+(const std::string& s) const;

        std::string path;
    };
}