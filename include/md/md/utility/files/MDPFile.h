#pragma once

#include <io/File.h>

namespace md {
    // Molecular dynamics parameter file
    struct MDPFile : public io::File {
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