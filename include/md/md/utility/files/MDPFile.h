#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace ausaxs::md {
    namespace detail {
        struct validate_mdp_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".mdp") {throw std::runtime_error("MDPFile::validate: File \"" + f->path() + "\" is not a molecular dynamics parameter file (.mdp).");}
            }
        };
    }

    // Molecular dynamics parameter file
    struct MDPFile : public io::detail::IValidatedFile<detail::validate_mdp_file> {
        using IValidatedFile::IValidatedFile;
    };
}