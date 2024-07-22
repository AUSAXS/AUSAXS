#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_gro_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".gro") {throw std::runtime_error("GROFile::validate: File \"" + f->path() + "\" is not a coordinate file (.gro).");}
            }
        };
    }

    // Coordinate file
    struct GROFile : public io::detail::IValidatedFile<detail::validate_gro_file> {
        using IValidatedFile::IValidatedFile;
        std::string get_unit_cell() const;
    };
}