#pragma once

#include <md/utility/files/File.h>

namespace md {
    // Trajectory file
    struct XTCFile : public detail::File {
        XTCFile() = default;
        XTCFile(const std::string& name) : File(name, "xtc") {}
        XTCFile(const char* name) : XTCFile(std::string(name)) {}
        ~XTCFile() override = default;
    };
}