#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_tpr_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".tpr") {throw std::runtime_error("TPRFile::validate: File \"" + f->path() + "\" is not a binary run input file (.tpr).");}
            }
        };
    }

    // Binary run input file
    struct TPRFile : public io::detail::IValidatedFile<detail::validate_tpr_file> {
        using IValidatedFile::IValidatedFile;
    };
}