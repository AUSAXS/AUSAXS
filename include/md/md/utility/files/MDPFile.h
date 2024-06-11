#pragma once

#include <md/utility/files/File.h>

namespace md {
    // Molecular dynamics parameter file
    struct MDPFile : public detail::File {
        MDPFile() = default;
        MDPFile(const std::string& name) : File(name, "mdp") {}
        MDPFile(const char* name) : MDPFile(std::string(name)) {}
        ~MDPFile() override = default;

        MDPFile& create() {
            File::create();
            return *this;
        }
    };
}