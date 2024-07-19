#pragma once

#include <io/File.h>

namespace md {
    // Checkpoint file
    struct CPTFile : public io::File {
        CPTFile() = default;
        CPTFile(const std::string& name) : File(name, "cpt") {}
        CPTFile(const char* name) : CPTFile(std::string(name)) {}
        ~CPTFile() override = default;
    };
}