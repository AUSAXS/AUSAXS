#pragma once

#include <md/utility/files/File.h>

namespace md {
    // Include topology file
    struct TXTFile : public detail::File {
        TXTFile() = default;
        TXTFile(const std::string& name) : File(name, "txt") {}
        TXTFile(const char* name) : TXTFile(std::string(name)) {}
        ~TXTFile() override = default;
    };
}