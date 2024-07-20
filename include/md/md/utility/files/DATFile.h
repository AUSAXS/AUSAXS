#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_dat_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".dat") {throw std::runtime_error("DATFile::validate: File \"" + f->path() + "\" is not a binary run input file (.dat).");}
            }
        };
    }

    // Binary run input file
    struct DATFile : public io::detail::IValidatedFile<detail::validate_dat_file> {
        using IValidatedFile::IValidatedFile;
    };
}