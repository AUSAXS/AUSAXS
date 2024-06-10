#pragma once

#include <md/utility/files/File.h>

namespace gmx {
    // Coordinate file
    struct GROFile : public detail::File {
        GROFile() = default;
        GROFile(const std::string& name) : File(name, "gro") {}
        GROFile(const char* name) : GROFile(std::string(name)) {}

        std::string get_unit_cell() const;
    };
}