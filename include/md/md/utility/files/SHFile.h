#pragma once

#include <io/File.h>

namespace md {
    // Binary run input file
    struct SHFile : public io::File {
        SHFile() = default;
        SHFile(const std::string& name) : File(name, "sh") {}
        SHFile(const char* name) : SHFile(std::string(name)) {}
        ~SHFile() override = default;
    };
}