#pragma once

#include <io/File.h>

namespace md {
    // Trajectory file
    struct XTCFile : public io::File {
        XTCFile() = default;
        XTCFile(const std::string& name) : File(name, "xtc") {}
        XTCFile(const char* name) : XTCFile(std::string(name)) {}
        ~XTCFile() override = default;
    };
}