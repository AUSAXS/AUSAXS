#pragma once

#include <md/utility/files/File.h>

namespace gmx {
    // Include topology file
    struct XVGFile : public detail::File {
        XVGFile() = default;
        XVGFile(const std::string& name) : File(name, "xvg") {}
        XVGFile(const char* name) : XVGFile(std::string(name)) {}
    };
}