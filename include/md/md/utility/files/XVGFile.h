#pragma once

#include <io/File.h>

namespace md {
    // Include topology file
    struct XVGFile : public io::File {
        XVGFile() = default;
        XVGFile(const std::string& name) : File(name, "xvg") {}
        XVGFile(const char* name) : XVGFile(std::string(name)) {}
        ~XVGFile() override = default;
    };
}