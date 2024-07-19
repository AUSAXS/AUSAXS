#pragma once

#include <io/File.h>

namespace md {
    // Include topology file
    struct TXTFile : public io::File {
        TXTFile() = default;
        TXTFile(const std::string& name) : File(name, "txt") {}
        TXTFile(const char* name) : TXTFile(std::string(name)) {}
        ~TXTFile() override = default;
    };
}