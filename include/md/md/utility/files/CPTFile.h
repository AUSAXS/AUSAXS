#pragma once

#include <md/utility/files/File.h>

namespace md {
    // Checkpoint file
    struct CPTFile : public detail::File {
        CPTFile() = default;
        CPTFile(const std::string& name) : File(name, "cpt") {}
        CPTFile(const char* name) : CPTFile(std::string(name)) {}
        ~CPTFile() override = default;
    };
}