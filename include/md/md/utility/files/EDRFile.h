#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace ausaxs::md {
    namespace detail {
        struct validate_edr_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".edr") {throw std::runtime_error("EDRFile::validate: File \"" + f->path() + "\" is not an energy file (.edr).");}
            }
        };
    }

    // Energy file
    struct EDRFile : public io::detail::IValidatedFile<detail::validate_edr_file> {
        using IValidatedFile::IValidatedFile;
    };
}