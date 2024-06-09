#pragma once

#include <utility/files/File.h>

namespace gmx {
    // Include topology file
    struct TXTFile : public detail::File {
        TXTFile() = default;
        TXTFile(const std::string& name) : File(name, "txt") {}
        TXTFile(const char* name) : TXTFile(std::string(name)) {}
    };
}