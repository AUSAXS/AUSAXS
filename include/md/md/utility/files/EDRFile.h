#pragma once

#include <io/File.h>

namespace md {
    // Energy file
    struct EDRFile : public io::File {
        EDRFile() = default;
        EDRFile(const std::string& name) : File(name, "edr") {}
        EDRFile(const char* name) : EDRFile(std::string(name)) {}
        ~EDRFile() override = default;
    };
}