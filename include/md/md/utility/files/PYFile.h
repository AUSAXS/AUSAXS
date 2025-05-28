#pragma once

#include <io/detail/IValidatedFile.h>
#include <utility/observer_ptr.h>

#include <stdexcept>

namespace ausaxs::md {
    namespace detail {
        struct validate_py_file {
            static void validate(observer_ptr<io::File> f) {
                if (f->extension() != ".py") {throw std::runtime_error("PYFile::validate: File \"" + f->path() + "\" is not a python file (.py).");}
            }
        };
    }

    // Python file
    struct PYFile : public io::detail::IValidatedFile<detail::validate_py_file> {
        using IValidatedFile::IValidatedFile;
    };
}