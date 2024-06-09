#pragma once

#include <utility/files/File.h>
#include <iomanip>

namespace gmx {
    // Molecular dynamics parameter file
    struct MDPFile : public detail::File {
        MDPFile() = default;
        MDPFile(const std::string& name) : File(name, "mdp") {}
        MDPFile(const char* name) : MDPFile(std::string(name)) {}

        MDPFile& create() {
            File::create();
            return *this;
        }
    };
}