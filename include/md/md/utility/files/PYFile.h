#pragma once

#include <io/File.h>

namespace md {
    // Binary run input file
    struct PYFile : public io::File {
        PYFile() = default;
        PYFile(const std::string& name) : File(name, "py") {}
        PYFile(const char* name) : PYFile(std::string(name)) {}
        ~PYFile() override = default;
    };
}