#pragma once

#include <utility/files/File.h>

namespace gmx {
    // Binary run input file
    struct SHFile : public detail::File {
        SHFile() = default;
        SHFile(const std::string& name) : File(name, "sh") {}
        SHFile(const char* name) : SHFile(std::string(name)) {}
    };
}