#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace ausaxs::md {
    namespace detail {
        struct validate_xtc_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".xtc") {throw std::runtime_error("XTCFile::validate: File \"" + f->path() + "\" is not a trajectory file (.xtc).");}
            }
        };
    }

    // Trajectory file
    struct XTCFile : public io::detail::IValidatedFile<detail::validate_xtc_file> {
        using IValidatedFile::IValidatedFile;
    };
}