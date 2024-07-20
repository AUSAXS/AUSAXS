#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_txt_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".txt") {throw std::runtime_error("TXTFile::validate: File \"" + f->path() + "\" is not a text file (.txt).");}
            }
        };
    }

    // Text file
    struct TXTFile : public io::detail::IValidatedFile<detail::validate_txt_file> {
        using IValidatedFile::IValidatedFile;
    };
}