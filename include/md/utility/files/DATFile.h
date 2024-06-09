#pragma once

#include <utility/files/File.h>

namespace gmx {
    // Binary run input file
    struct DATFile : public detail::File {
        DATFile() = default;
        DATFile(const std::string& name) : File(name, "dat") {}
        DATFile(const char* name) : DATFile(std::string(name)) {}
    };
}