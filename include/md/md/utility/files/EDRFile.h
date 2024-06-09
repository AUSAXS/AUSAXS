#pragma once

#include <utility/files/File.h>

namespace gmx {
    // Energy file
    struct EDRFile : public detail::File {
        EDRFile() = default;
        EDRFile(const std::string& name) : File(name, "edr") {}
        EDRFile(const char* name) : EDRFile(std::string(name)) {}
    };
}