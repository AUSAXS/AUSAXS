#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_cpt_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".cpt") {throw std::runtime_error("CPTFile::validate: File \"" + f->path() + "\" is not a checkpoint file (.cpt).");}
            }
        };
    }

    // Checkpoint file
    struct CPTFile : public io::detail::IValidatedFile<detail::validate_cpt_file> {
        using IValidatedFile::IValidatedFile;
    };
}