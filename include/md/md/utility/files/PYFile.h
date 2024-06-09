#pragma once

#include <utility/files/File.h>

namespace gmx {
    // Binary run input file
    struct PYFile : public detail::File {
        PYFile() = default;
        PYFile(const std::string& name) : File(name, "py") {}
        PYFile(const char* name) : PYFile(std::string(name)) {}
    };
}