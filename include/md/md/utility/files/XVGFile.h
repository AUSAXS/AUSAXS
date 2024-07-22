#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace md {
    namespace detail {
        struct validate_xvg_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".xvg") {throw std::runtime_error("XVGFile::validate: File \"" + f->path() + "\" is not a GROMACS data output file (.xvg).");}
            }
        };
    }

    // Molecular dynamics parameter file
    struct XVGFile : public io::detail::IValidatedFile<detail::validate_xvg_file> {
        using IValidatedFile::IValidatedFile;
    };
}