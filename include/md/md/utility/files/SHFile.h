#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace ausaxs::md {
    namespace detail {
        struct validate_sh_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".sh") {throw std::runtime_error("SHFile::validate: File \"" + f->path() + "\" is not a shell script (.sh).");}
            }
        };
    }

    // Shell script file
    struct SHFile : public io::detail::IValidatedFile<detail::validate_sh_file> {
        using IValidatedFile::IValidatedFile;
    };
}